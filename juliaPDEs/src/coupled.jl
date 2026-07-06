# ── CoupledPDESystem — M linearly-coupled fields, one shared grid ────────────
#
#   ∂uₖ/∂t = αₖ Δuₖ − Σⱼ Kₖⱼ(uₖ − uⱼ)      k = 1, …, M
#
# on the Cartesian box ∏ᵢ [aᵢ, bᵢ] in ℝᴺ, homogeneous Dirichlet (every field
# = 0) on the boundary. `K` is an M×M linear exchange-rate matrix between
# fields — e.g. M plates at different temperatures, each diffusing internally
# while conducting heat into whichever other plates it's coupled to (`K[i][j]`
# nonzero); `K` all-zero recovers M independent, uncoupled fields solved
# together purely for the convenience of one shared `Ind`-indexed system.
#
# This generalises the earlier 2-field-only `CoupledHeatEquation` (still kept
# below as a convenience constructor) to any number M of named fields,
# stacked into ONE sparse linear system per time step instead of time-stepping
# each field with its own independent solve — the only way to represent a
# term that appears in more than one equation at the same implicit time
# level. The M interior fields are concatenated into one vector
# w = [u₁ⁿ⁺¹; u₂ⁿ⁺¹; …; u_Mⁿ⁺¹] (via the same flat layout `Ind` uses) and a
# single (Mn)×(Mn) block system is assembled per step:
#
#   Aᵢᵢ = I + dt(αᵢ·Lap + (Σⱼ Kᵢⱼ)·I)          diagonal block
#   Aᵢⱼ = −dt·Kᵢⱼ·I                             off-diagonal block  (i ≠ j)
#
# where `Lap = nd_laplacian(...)` is the same sparse operator (Lap ≈ −Δ) that
# `PoissonEquation` (poisson.jl) and `HeatEquation.solve_implicit` (heat.jl)
# already build — backward Euler, so the scheme stays unconditionally stable
# regardless of K. The block matrix is factorised once via `lu` and reused
# for every time step.
#
# Naming (`vars::NTuple{M,Symbol}`, the flat `Ind`-based layout) mirrors the
# indexing scheme from the exploratory generic-PDE project this multi-field
# support is modelled on (see `Ind` in types.jl).
#
Base.@kwdef struct CoupledPDESystem{M,N} <: ParabolicProblem
  grid::Grid{N}
  vars::NTuple{M,Symbol}                                             # field names, e.g. (:u, :v, :w)
  α::NTuple{M,Float64}                                               # diffusivity of each field
  # NOTE: the default below references `vars` (a previously-listed field), not the type
  # parameter `M` — `@kwdef` default expressions are evaluated in the outer constructor
  # before M is resolved, so `M` itself is not in scope here.
  K::NTuple{M,NTuple{M,Float64}} = ntuple(_ -> ntuple(_ -> 0.0, length(vars)), length(vars)) # coupling matrix, default = no coupling
  f_init::NTuple{M,Function}                                         # initial condition, one closure per field
  Nt::Int = 1000                                                     # number of time steps (dt = T / Nt)
  T::Float64 = 1.0                                                   # final simulation time
end

function solve(p::CoupledPDESystem{M,N}) where {M,N}
  space = p.grid.space
  npts  = p.grid.numgrid

  d         = ntuple(i -> (space[i][end] - space[i][1]) / (npts[i] - 1), N)
  inner_n   = ntuple(i -> npts[i] - 2, N)
  inner_idx = ntuple(i -> 2:(npts[i]-1), N)
  all(>(0), inner_n) || error("CoupledPDESystem needs ≥3 grid points per axis (got npts=$npts).")
  n = prod(inner_n)

  dt = p.T / p.Nt

  # 1. Shared spatial operator and one-time (Mn)×(Mn) block-system factorisation.
  Lap   = nd_laplacian(d, inner_n)
  Iden  = sparse(I, n, n)
  koff  = ntuple(i -> sum(p.K[i][j] for j in 1:M if j != i), M)
  A = vcat((hcat((i == j ? Iden + dt * (p.α[i] * Lap + koff[i] * Iden) : -dt * p.K[i][j] * Iden
                  for j in 1:M)...) for i in 1:M)...)
  Afact = lu(A)

  # 2. Initial condition per field, then stamp homogeneous Dirichlet (= 0) on
  #    every boundary face — f_init is only trusted on the interior.
  fields = ntuple(_ -> zeros(Float64, npts), M)
  for k in 1:M
    for I in CartesianIndices(fields[k])
      coords = ntuple(dim -> space[dim][I[dim]], N)
      fields[k][I] = p.f_init[k](coords...)
    end
    for i in 1:N
      selectdim(fields[k], i, 1)       .= 0.0
      selectdim(fields[k], i, npts[i]) .= 0.0
    end
  end

  ints = [vec(fields[k][inner_idx...]) for k in 1:M]

  # 3. Time stepping — one sparse solve per step advances ALL M fields together.
  for _ in 1:p.Nt
    w = Afact \ vcat(ints...)
    for k in 1:M
      ints[k] = w[(k-1)*n+1:k*n]
    end
  end

  for k in 1:M
    fields[k][inner_idx...] = reshape(ints[k], inner_n)
  end

  return CoupledPDESolution(space, p.vars, fields, p.T, p)
end

# ── CoupledHeatEquation — the original 2-field case ──────────────────────────
#
#   ∂u/∂t = α₁ Δu − κ(u − v)
#   ∂v/∂t = α₂ Δv − κ(v − u)
#
# Kept as a convenience constructor: it now builds and delegates to a
# `CoupledPDESystem{2,N}` with `vars = (:u, :v)` and `K = [0 κ; κ 0]`, so the
# 2-field special case is a literal instance of the general M-field system
# rather than a hand-maintained parallel implementation — `sol.u`/`sol.v`
# keep working exactly as before via `CoupledPDESolution`'s `getproperty`.
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
  return solve(CoupledPDESystem(
    grid   = p.grid,
    vars   = (:u, :v),
    α      = p.α,
    K      = ((0.0, p.κ), (p.κ, 0.0)),
    f_init = p.f_init,
    Nt     = p.Nt,
    T      = p.T,
  ))
end
