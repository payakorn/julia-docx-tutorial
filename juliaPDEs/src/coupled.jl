# ── Coupled Heat Equations — two linearly-coupled diffusing fields ───────────
#
#   ∂u/∂t = α₁ Δu − κ(u − v)
#   ∂v/∂t = α₂ Δv − κ(v − u)
#
# on the Cartesian box ∏ᵢ [aᵢ, bᵢ] in ℝᴺ, homogeneous Dirichlet (u=v=0) on the
# boundary. κ is a linear exchange rate between the two fields — e.g. two
# plates at different temperatures, each diffusing internally while also
# conducting heat into one another.
#
# This is the smallest example of solving a MULTI-FIELD system as ONE sparse
# linear system per time step, instead of time-stepping each field with its
# own independent solve (which cannot represent a term that appears in both
# equations at the same implicit time level). Stack the two interior fields
# into one vector w = [uⁿ⁺¹; vⁿ⁺¹] and assemble a single (2n)×(2n) block
# system per step:
#
#   [ I + dt(α₁·Lap + κI)      −dt·κ·I          ] [uⁿ⁺¹]   [uⁿ]
#   [      −dt·κ·I          I + dt(α₂·Lap + κI) ] [vⁿ⁺¹] = [vⁿ]
#
# where `Lap = nd_laplacian(...)` is the same sparse operator (Lap ≈ −Δ) that
# `PoissonEquation` (poisson.jl) and `HeatEquation.solve_implicit` (heat.jl)
# already build — backward Euler, so the scheme stays unconditionally stable
# regardless of κ. The block matrix is factorised once via `lu` and reused
# for every time step.
#
Base.@kwdef struct CoupledHeatEquation{N,F1,F2} <: ParabolicProblem
  grid::Grid{N}
  α::NTuple{2,Float64} = (0.01, 0.02)   # diffusivity of field 1, field 2
  κ::Float64           = 0.5            # linear exchange rate between fields
  f_init::Tuple{F1,F2}                  # initial condition, one closure per field
  Nt::Int              = 1000           # number of time steps (dt = T / Nt)
  T::Float64           = 1.0            # final simulation time
end

function solve(p::CoupledHeatEquation{N,F1,F2}) where {N,F1,F2}
  space = p.grid.space
  npts  = p.grid.numgrid

  d         = ntuple(i -> (space[i][end] - space[i][1]) / (npts[i] - 1), N)
  inner_n   = ntuple(i -> npts[i] - 2, N)
  inner_idx = ntuple(i -> 2:(npts[i]-1), N)
  all(>(0), inner_n) || error("CoupledHeatEquation needs ≥3 grid points per axis (got npts=$npts).")
  n = prod(inner_n)

  dt = p.T / p.Nt

  # 1. Shared spatial operator and one-time block-system factorisation.
  Lap  = nd_laplacian(d, inner_n)
  Iden = sparse(I, n, n)
  A11  = Iden + dt * (p.α[1] * Lap + p.κ * Iden)
  A22  = Iden + dt * (p.α[2] * Lap + p.κ * Iden)
  Aoff = -dt * p.κ * Iden
  Afact = lu([A11 Aoff; Aoff A22])

  # 2. Initial condition, then stamp homogeneous Dirichlet (u=v=0) on every
  #    boundary face — f_init is only trusted on the interior.
  u = zeros(Float64, npts)
  v = zeros(Float64, npts)
  for I in CartesianIndices(u)
    coords = ntuple(dim -> space[dim][I[dim]], N)
    u[I] = p.f_init[1](coords...)
    v[I] = p.f_init[2](coords...)
  end
  for i in 1:N
    selectdim(u, i, 1)       .= 0.0
    selectdim(u, i, npts[i]) .= 0.0
    selectdim(v, i, 1)       .= 0.0
    selectdim(v, i, npts[i]) .= 0.0
  end

  u_int = vec(u[inner_idx...])
  v_int = vec(v[inner_idx...])

  # 3. Time stepping — one sparse solve per step advances BOTH fields together.
  for _ in 1:p.Nt
    w = Afact \ vcat(u_int, v_int)
    u_int = w[1:n]
    v_int = w[n+1:end]
  end

  u[inner_idx...] = reshape(u_int, inner_n)
  v[inner_idx...] = reshape(v_int, inner_n)

  return CoupledPDESolution(space, u, v, p.T, p)
end
