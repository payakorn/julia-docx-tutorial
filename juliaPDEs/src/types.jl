# ── Abstract type hierarchy ───────────────────────────────────────────────────
#
#   PDEProblem
#   ├── ParabolicProblem    (time-dependent, diffusion-dominated — e.g. Heat)
#   ├── HyperbolicProblem   (time-dependent, wave-like — e.g. Wave)
#   ├── EllipticProblem     (steady-state — e.g. Poisson)
#   └── IncompressibleNSProblem  (Navier-Stokes)
#
abstract type PDEProblem end
abstract type ParabolicProblem <: PDEProblem end
abstract type HyperbolicProblem <: PDEProblem end
abstract type EllipticProblem <: PDEProblem end
abstract type IncompressibleNSProblem <: PDEProblem end

# ── PDESolution — works for any number of dimensions ─────────────────────────
#
#   Subtypes AbstractArray{T,N} so the full Julia array interface comes
#   automatically: length, size, iterate, broadcasting, slicing, sum, …
#
#   grid  — one coordinate vector per spatial (or time) dimension
#   u     — the field values, same shape as the Cartesian product of grid
#   t     — final simulation time
#   problem — the original problem struct (preserves parameters for plotting)
#
struct PDESolution{T,N} <: AbstractArray{T,N}
  grid::NTuple{N,Vector{Float64}}
  u::Array{T,N}
  t::Float64
  problem::PDEProblem
end


Base.@kwdef struct Grid{N}
  a::NTuple{N,Float64}
  b::NTuple{N,Float64}
  stepsize::NTuple{N,Float64}

  # 1. Use broadcasting with Tuple to ensure the type matches NTuple{N, Int}
  numgrid::NTuple{N,Int} = Int.(ceil.((b .- a) ./ stepsize))

  # 2. Use ntuple() for clean, type-stable generation of the spaces
  space::NTuple{N,Vector{Float64}} = ntuple(i -> collect(range(a[i], b[i], length=numgrid[i])), length(a))
end


struct TestGrid{N,F}
  grid::Grid{N}
  bc1::NTuple{N,Function}
  bc2::NTuple{N,Function}
  f_init::F
end

# ── Ind — flat linear indexing across M named fields sharing one grid ────────
#
#   Mirrors the indexing scheme from the exploratory generic-PDE project this
#   package's multi-field support is modelled on: rather than hardcoding
#   separate `u`, `v`, … fields on a solution/problem struct, stack the M
#   interior fields into ONE vector (as `CoupledPDESystem.solve` does for its
#   block linear system) and use `vars::NTuple{M,Symbol}` to look positions
#   up by name.
#
#   `Ind(I::Int)`        → (which variable, CartesianIndex within that field)
#   `Ind(var, I...)`     → flat linear index of `var[I...]` in the stacked vector
#
#   Unlike the exploratory version this package draws on — which stored one
#   `Grid` per field — every field here shares the single `Grid{N}` of the
#   `CoupledPDESystem` it belongs to, since a coupled system lives on one
#   common domain. `NumGrid` still carries one shape per field (all equal)
#   so the indexing arithmetic reads the same as the source it mirrors.
#
struct Ind{M,N}
  vars::NTuple{M,Symbol}
  NumGrid::NTuple{M,NTuple{N,Int}}
end

Ind(vars::NTuple{M,Symbol}, numgrid::NTuple{N,Int}) where {M,N} =
  Ind(vars, ntuple(_ -> numgrid, M))

function Base.getindex(d::Ind{M,N}, I::Int) where {M,N}
  n = prod(d.NumGrid[1])
  nvar = I % n == 0 ? I ÷ n : I ÷ n + 1
  ngrid = I - (nvar - 1) * n
  return d.vars[nvar], CartesianIndices(d.NumGrid[1])[ngrid]
end

function Base.getindex(d::Ind{M,N}, var::Symbol, I...) where {M,N}
  length(I) == N || error("Ind: dimension mismatch (expected $N indices, got $(length(I)))")
  nvar = findfirst(==(var), d.vars)
  nvar === nothing && error("Ind: unknown variable $var (known: $(d.vars))")
  return prod(d.NumGrid[1]) * (nvar - 1) + LinearIndices(d.NumGrid[1])[I...]
end

# ── CoupledPDESolution — M fields sharing one grid ────────────────────────────
#
#   A second, minimal solution type rather than reusing `PDESolution{T,N}`:
#   that type subtypes `AbstractArray{T,N}` around a single homogeneous array,
#   so bolting extra fields onto it would change what indexing/iterating a
#   `PDESolution` means for every existing solver. Coupled problems (M fields
#   solved together as one linear system, e.g. `CoupledPDESystem`) get their
#   own lightweight struct instead.
#
#   `vars`/`fields` generalise the earlier fixed `u`/`v` design to any number
#   of named fields; `sol.u`, `sol.v`, … keep working via `getproperty` for
#   any name present in `vars` (the `CoupledHeatEquation` 2-field case looks
#   identical to before — it now just happens to be `M == 2`).
#
struct CoupledPDESolution{M,N}
  grid::NTuple{N,Vector{Float64}}
  vars::NTuple{M,Symbol}
  fields::NTuple{M,Array{Float64,N}}
  t::Float64
  problem::PDEProblem
end

function Base.getproperty(s::CoupledPDESolution, sym::Symbol)
  vars = getfield(s, :vars)
  i = findfirst(==(sym), vars)
  i !== nothing && return getfield(s, :fields)[i]
  g = getfield(s, :grid)
  sym === :x && return g[1]
  if sym === :y
    length(g) >= 2 || error("$(typeof(s)) has no y dimension")
    return g[2]
  end
  if sym === :z
    length(g) >= 3 || error("$(typeof(s)) has no z dimension")
    return g[3]
  end
  return getfield(s, sym)
end

function Base.show(io::IO, ::MIME"text/plain", s::CoupledPDESolution{M,N}) where {M,N}
  dims = join(["$(length(g)) pts" for g in s.grid], " × ")
  println(io, "CoupledPDESolution{$(M) fields, $(N)D}")
  println(io, "  problem : $(typeof(s.problem))")
  println(io, "  size    : $(size(s.fields[1])) (× $(M) fields: $(join(s.vars, ", ")))")
  println(io, "  t       : $(s.t)")
  print(io, "  grid    : $(dims)")
end

# ── AbstractArray interface ───────────────────────────────────────────────────
Base.size(s::PDESolution) = size(s.u)
Base.getindex(s::PDESolution, I...) = s.u[I...]
Base.setindex!(s::PDESolution, v, I...) = (s.u[I...] = v)
Base.IndexStyle(::Type{<:PDESolution}) = IndexCartesian()

# ── Convenience properties: sol.x, sol.y, sol.z ──────────────────────────────
function Base.getproperty(s::PDESolution, sym::Symbol)
  g = getfield(s, :grid)
  sym === :x && return g[1]
  if sym === :y
    length(g) >= 2 || error("$(typeof(s)) has no y dimension")
    return g[2]
  end
  if sym === :z
    length(g) >= 3 || error("$(typeof(s)) has no z dimension")
    return g[3]
  end
  return getfield(s, sym)
end

function Base.propertynames(s::PDESolution{T,N}, private::Bool=false) where {T,N}
  names = (:grid, :u, :t, :problem, :x)
  N >= 2 && (names = (names..., :y))
  N >= 3 && (names = (names..., :z))
  return names
end

# ── Pretty printing — suppress the default "print every element" dump ─────────
function Base.show(io::IO, ::MIME"text/plain", s::PDESolution{T,N}) where {T,N}
  dims = join(["$(length(g)) pts" for g in s.grid], " × ")
  println(io, "PDESolution{$(N)D, $(T)}")
  println(io, "  problem : $(typeof(s.problem))")
  println(io, "  size    : $(size(s.u))")
  println(io, "  t       : $(s.t)")
  print(io, "  grid    : $(dims)")
end

Base.show(io::IO, s::PDESolution{T,N}) where {T,N} =
  print(io, "PDESolution{$(N)D}($(join(size(s.u), "×")), t=$(s.t))")


# ── Grid builders ─────────────────────────────────────────────────────────────
#
# Two conventions show up across the PDE solvers in this package:
#
#  • interior_grid — for Dirichlet problems whose boundary values are baked in
#    and not stored. We solve only for the N interior unknowns; x[1] sits at
#    Δx, not 0. Used by Heat and Poisson.
#
#  • endpoint_grid — for problems whose boundary values are part of the array
#    (e.g. wave equation, where the leap-frog update reads u[1] and u[end]).
#    The grid spans the closed interval [0, L] with N evenly-spaced points.
#
# Both helpers return (Δx, x). Centralising the formulae here means a fix in
# one place propagates to every solver.

"""
    interior_grid(a, b, N) -> (dx, x)

Uniform grid of `N` interior points on the open interval `(a, b)`.
`dx = (b - a) / (N + 1)` and `x[i] = a + i·dx`. The boundary nodes
`a` and `b` are intentionally not stored.

    interior_grid(L, N) = interior_grid(0, L, N)
"""
function interior_grid(a::Real, b::Real, N::Integer)
  dx = (b - a) / (N + 1)
  x = collect(range(a + dx, b - dx, length=N))
  return dx, x
end
interior_grid(L::Real, N::Integer) = interior_grid(zero(L), L, N)

"""
    endpoint_grid(a, b, N) -> (dx, x)

Uniform grid of `N` points on the closed interval `[a, b]`, endpoints
included. `dx = (b - a) / (N - 1)` and `x = [a, a+dx, …, b]`.

    endpoint_grid(L, N) = endpoint_grid(0, L, N)
"""
function endpoint_grid(a::Real, b::Real, N::Integer)
  dx = (b - a) / (N - 1)
  x = collect(range(a, b, length=N))
  return dx, x
end
endpoint_grid(L::Real, N::Integer) = endpoint_grid(zero(L), L, N)
