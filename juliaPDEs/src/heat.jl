# ── Heat equation — dimension-free, TestGrid-driven ──────────────────────────
#
#   ∂u/∂t = α Δu          on the Cartesian box  ∏ᵢ [aᵢ, bᵢ]    in ℝᴺ
#
#   Dirichlet boundary values are supplied by the TestGrid:
#       lower face of axis i :  u = bc1[i](t, x_other...)
#       upper face of axis i :  u = bc2[i](t, x_other...)
#   where `x_other...` are the coordinates on the (N−1) axes orthogonal to i,
#   in their natural order. For N = 1 the call degenerates to bc?[1](t).
#
#   u(x, 0) = f_init(x...)        (TestGrid.f_init)
#
# Parameter `N` is the spatial dimension (inherited from the TestGrid);
# `F` is the initial-condition closure type, so Julia specialises the solver
# per f_init the user passes through TestGrid.
#
Base.@kwdef struct HeatEquation{N,F} <: ParabolicProblem
  testgrid::TestGrid{N,F}
  Nt::Int = 1000        # number of time steps (dt = T / Nt)
  T::Float64 = 1.0         # final simulation time
  α::Float64 = 0.01        # diffusivity
end

# Overwrite every Dirichlet face of `u` with the values returned by bc1/bc2 at
# time `t`. For each axis i, walk the (N−1)-D face slice and evaluate the BC
# function on the orthogonal coordinates.
@inline function _apply_dirichlet_bcs!(u::AbstractArray{Float64,N},
  bc1::NTuple{N,Function},
  bc2::NTuple{N,Function},
  space::NTuple{N,Vector{Float64}},
  t::Float64) where {N}
  for i in 1:N
    face_lo = selectdim(u, i, 1)
    face_hi = selectdim(u, i, size(u, i))
    for J in CartesianIndices(face_lo)
      # k iterates 1..N−1 over the face axes; map back to the original axis
      # by skipping `i`.
      coords_other = ntuple(k -> space[k < i ? k : k + 1][J[k]], N - 1)
      face_lo[J] = bc1[i](t, coords_other...)
      face_hi[J] = bc2[i](t, coords_other...)
    end
  end
  return u
end

function solve(p::HeatEquation{N,F}) where {N,F}
  tg = p.testgrid
  space = tg.grid.space
  npts = tg.grid.numgrid

  # 1. Actual per-axis spacing from the closed-interval grid. Computed from
  #    `space` rather than `tg.grid.stepsize` because `numgrid` rounds up and
  #    the requested stepsize need not evenly divide (b − a).
  d = ntuple(i -> (space[i][end] - space[i][1]) / (npts[i] - 1), N)

  # 2. Forward-Euler stability: Σᵢ α dt / dxᵢ² ≤ 1/2.
  dt = p.T / p.Nt
  r = ntuple(i -> p.α * dt / d[i]^2, N)
  sum(r) ≤ 0.5 || @warn "HeatEquation stability: Σrᵢ = $(sum(r)) > 0.5 — increase Nt."

  # 3. Allocate on the endpoint grid and seed from f_init.
  u = zeros(Float64, npts)
  for I in CartesianIndices(u)
    coords = ntuple(dim -> space[dim][I[dim]], N)
    u[I] = tg.f_init(coords...)
  end

  # 4. Stamp the initial Dirichlet boundary values at t = 0 (they win over
  #    whatever f_init produced on the faces).
  _apply_dirichlet_bcs!(u, tg.bc1, tg.bc2, space, 0.0)

  u_new = copy(u)

  inner_range = CartesianIndices(ntuple(i -> 2:(npts[i]-1), N))
  e_off = ntuple(i -> CartesianIndex(ntuple(j -> j == i ? 1 : 0, N)), N)

  # 5. Time-stepping loop — forward-Euler + (2N+1)-point Laplacian. Boundary
  #    values are refreshed at the new time level after each interior sweep.
  for n in 1:p.Nt
    for I in inner_range
      laplacian = 0.0
      for i in 1:N
        @inbounds laplacian += r[i] * (u[I+e_off[i]] - 2.0 * u[I] + u[I-e_off[i]])
      end
      @inbounds u_new[I] = u[I] + laplacian
    end

    _apply_dirichlet_bcs!(u_new, tg.bc1, tg.bc2, space, n * dt)

    u, u_new = u_new, u
  end

  return PDESolution(space, u, p.T, p)
end

# ── Implicit time stepping: θ-method (Crank-Nicolson by default) ─────────────
#
#   (I − α dt θ L) u^{n+1} = (I + α dt (1−θ) L) u^n
#                            + α dt (θ g^{n+1} + (1−θ) g^n)
#
# where L is the sparse Laplacian on the interior (built once via the
# Kronecker-sum helper in poisson.jl, sign-flipped because `nd_laplacian`
# returns −Δ) and `g` is the boundary lift vector — known Dirichlet values
# at the faces, scaled by 1/dxᵢ², that would otherwise sit on the RHS of
# the central-difference stencil for nodes adjacent to a wall.
#
# θ = 1/2 is Crank-Nicolson (2nd order in time, unconditionally stable).
# θ = 1.0 is backward Euler (1st order, also unconditionally stable, more
# dissipative — useful when you need strong damping of high-frequency modes).
#
function solve_implicit(p::HeatEquation{N,F}; θ::Float64=0.5) where {N,F}
  tg = p.testgrid
  space = tg.grid.space
  npts = tg.grid.numgrid

  d = ntuple(i -> (space[i][end] - space[i][1]) / (npts[i] - 1), N)
  inner_n = ntuple(i -> npts[i] - 2, N)
  inner_idx_full = ntuple(i -> 2:(npts[i]-1), N)

  all(>(0), inner_n) || error("solve_implicit needs ≥3 grid points per axis (got npts=$npts).")

  dt = p.T / p.Nt

  # 1. Assemble the operator. nd_laplacian returns −Δ; flip the sign for +Δ.
  L = -nd_laplacian(d, inner_n)
  Imat = sparse(I, prod(inner_n), prod(inner_n))
  A = Imat - (p.α * dt * θ) * L
  B = Imat + (p.α * dt * (1 - θ)) * L
  Afact = lu(A)                          # one factorisation, reused every step

  # 2. Initialise the full field and stamp BCs at t = 0.
  u = zeros(Float64, npts)
  for I in CartesianIndices(u)
    coords = ntuple(dim -> space[dim][I[dim]], N)
    u[I] = tg.f_init(coords...)
  end
  _apply_dirichlet_bcs!(u, tg.bc1, tg.bc2, space, 0.0)

  # 3. Interior solution vector (linear-algebra view).
  u_int = vec(u[inner_idx_full...])
  g_n = _boundary_lift(tg, space, d, inner_n, 0.0)

  # 4. Time stepping.
  for n in 1:p.Nt
    t_np1 = n * dt
    g_np1 = _boundary_lift(tg, space, d, inner_n, t_np1)
    rhs = B * u_int .+ (p.α * dt) .* (θ .* g_np1 .+ (1 - θ) .* g_n)
    u_int = Afact \ rhs
    g_n = g_np1
  end

  # 5. Pack the interior solution back and refresh boundaries at T.
  u[inner_idx_full...] = reshape(u_int, inner_n)
  _apply_dirichlet_bcs!(u, tg.bc1, tg.bc2, space, p.T)

  return PDESolution(space, u, p.T, p)
end

# Vector g such that  Δu_interior = L · u_interior + g  encodes the known
# Dirichlet boundary values at time `t` lifted onto interior nodes adjacent
# to each face. For axis i, the interior node at face-index 1 sees the lower
# boundary value, contributing bc1[i](t, x_other...) / dxᵢ²; symmetrically
# for the upper face via bc2[i].
function _boundary_lift(tg::TestGrid{N,F},
  space::NTuple{N,Vector{Float64}},
  d::NTuple{N,Float64},
  inner_n::NTuple{N,Int},
  t::Float64) where {N,F}
  g = zeros(Float64, inner_n)
  for i in 1:N
    face_lo = selectdim(g, i, 1)
    face_hi = selectdim(g, i, inner_n[i])
    for J in CartesianIndices(face_lo)
      # `J` is in INTERIOR coords on the (N−1) face axes. Map face axis k
      # back to the original axis (skipping i) and shift by +1 to land in
      # the FULL-grid coordinate vector (interior begins at full index 2).
      coords_other = ntuple(k -> space[k < i ? k : k + 1][J[k]+1], N - 1)
      face_lo[J] += tg.bc1[i](t, coords_other...) / d[i]^2
      face_hi[J] += tg.bc2[i](t, coords_other...) / d[i]^2
    end
  end
  return vec(g)
end
