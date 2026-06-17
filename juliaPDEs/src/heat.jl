# ── Heat equation — dimension-free ────────────────────────────────────────────
#
#   ∂u/∂t = α Δu          on the Cartesian box  ∏ᵢ [aᵢ, bᵢ]    in ℝᴺ
#   u = 0  on ∂Ω
#   u(x, 0) = f_init(x...)
#
# Parameter `N` is the spatial dimension (inferred from N_grid). Parameter `F`
# is the type of the initial-condition closure, so Julia specialises the solver
# for each f_init the user passes.
#
Base.@kwdef struct HeatEquation{N,F} <: ParabolicProblem
  N_grid::NTuple{N,Int} = (200,)                              # interior points per axis (defaults to 1D)
  a::NTuple{N,Float64} = ntuple(_ -> 0.0, length(N_grid))     # lower bound per axis
  b::NTuple{N,Float64} = ntuple(_ -> 1.0, length(N_grid))     # upper bound per axis
  Nt::Int = 1000                                              # number of time steps (dt = T / Nt)
  T::Float64 = 1.0                                            # final simulation time
  α::Float64 = 0.01                                           # diffusivity
  f_init::F                                                   # required — no default (was: x -> exp(-100*(x-0.5)^2))
end

function solve(p::HeatEquation{N,F}; save_every::Integer=0,
               save_dir::Union{Nothing,AbstractString}=nothing) where {N,F}
  # 1. Build the per-axis interior grids via the shared helper.
  grids = ntuple(i -> interior_grid(p.a[i], p.b[i], p.N_grid[i]), N)
  d = ntuple(i -> grids[i][1], N)
  axes_coords = ntuple(i -> grids[i][2], N)

  # 2. Time step from user-supplied Nt, with a stability check.
  dt = p.T / p.Nt
  r = ntuple(i -> p.α * dt / d[i]^2, N)
  sum(r) ≤ 0.5 || @warn "HeatEquation stability: Σrᵢ = $(sum(r)) > 0.5 — increase Nt."
  nsteps = p.Nt

  # 3. Initialize the field by sampling f_init on the Cartesian product grid.
  u = zeros(Float64, p.N_grid)
  for I in CartesianIndices(u)
    coords = ntuple(dim -> axes_coords[dim][I[dim]], N)
    u[I] = p.f_init(coords...)
  end

  # 4. Enforce Dirichlet zero boundaries ONCE at initialization.
  for i in 1:N
    selectdim(u, i, 1) .= 0.0
    selectdim(u, i, p.N_grid[i]) .= 0.0
  end

  u_new = copy(u)   # Clone to preserve initial boundaries

  # Setup inner range and stencil offsets — both fixed at compile time for known N.
  inner_range = CartesianIndices(ntuple(i -> 2:(p.N_grid[i]-1), N))
  e = ntuple(i -> CartesianIndex(ntuple(j -> j == i ? 1 : 0, N)), N)

  # Optional history writer — created only when save_every > 0 (fast path otherwise).
  writer = nothing
  if save_every > 0
    dir    = something(save_dir, joinpath("runs", default_run_name(p)))
    writer = SolutionWriter(dir; problem=p, grid=axes_coords, dt=dt, shape=p.N_grid)
    save_step!(writer, u, 0, 0, 0.0)          # initial condition = frame 0
  end

  # 5. Main time-stepping loop — forward-Euler + (2N+1)-point Laplacian.
  frame = 0
  for step in 1:nsteps
    for I in inner_range
      laplacian = 0.0
      for i in 1:N
        @inbounds laplacian += r[i] * (u[I+e[i]] - 2.0 * u[I] + u[I-e[i]])
      end
      @inbounds u_new[I] = u[I] + laplacian
    end
    # O(1) pointer swap instead of copying data (u .= u_new)
    u, u_new = u_new, u

    if writer !== nothing && step % save_every == 0
      frame += 1
      save_step!(writer, u, frame, step, step * dt)
    end
  end

  # Capture the final state (if not already on a save boundary), then write meta once.
  if writer !== nothing
    if nsteps % save_every != 0
      frame += 1
      save_step!(writer, u, frame, nsteps, nsteps * dt)
    end
    write_meta!(writer; params = _problem_params(p))
  end

  return PDESolution(axes_coords, u, p.T, p)
end
