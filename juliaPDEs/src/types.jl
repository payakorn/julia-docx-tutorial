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


# ── Types and Structs ─────────────────────────────────────────────────────────

# abstract type ParabolicProblem end

# Fully parameterized solution struct to guarantee type stability
# struct PDESolution{A, U, P}
#     axes::A
#     u::U
#     T::Float64
#     prob::P
# end

# Parameterize the function type `F` to avoid dynamic dispatch overhead
Base.@kwdef struct HeatEquationND{N,F} <: ParabolicProblem
    N_grid::NTuple{N,Int}
    Nt::Int = 1000                          # number of time steps (dt = T / Nt)
    α::Float64 = 0.01
    L::NTuple{N,Float64} = ntuple(_ -> 1.0, length(N_grid))
    T::Float64 = 1.0
    f_init::F
end


# ── Generic N-Dimensional Solver ──────────────────────────────────────────────

function solve(p::HeatEquationND{N,F}) where {N,F}
    # 1. Calculate grid spacings
    d = ntuple(i -> p.L[i] / (p.N_grid[i] + 1), N)

    # Generate physical coordinates for each axis
    axes_coords = ntuple(i -> collect(range(d[i], p.L[i] - d[i], length=p.N_grid[i])), N)

    # 2. Time step from user-supplied Nt, with a stability check
    dt = p.T / p.Nt
    r  = ntuple(i -> p.α * dt / d[i]^2, N)
    sum(r) ≤ 0.5 || @warn "HeatEquationND stability: Σrᵢ = $(sum(r)) > 0.5 — increase Nt."
    nsteps = p.Nt

    # 3. Initialize grid and map indices to physical coordinates
    u = zeros(Float64, p.N_grid)
    for I in CartesianIndices(u)
        coords = ntuple(dim -> axes_coords[dim][I[dim]], N)
        u[I] = p.f_init(coords...)
    end

    # 4. Enforce Dirichlet zero boundaries ONCE at initialization
    for i in 1:N
        selectdim(u, i, 1) .= 0.0
        selectdim(u, i, p.N_grid[i]) .= 0.0
    end

    u_new = copy(u) # Clone to preserve initial boundaries

    # Setup inner range to safely skip boundaries during the loop
    inner_range = CartesianIndices(ntuple(i -> 2:(p.N_grid[i]-1), N))

    # Pre-compute Cartesian offsets for the stencil (e.g., ±1 in X, Y, Z)
    e = ntuple(i -> CartesianIndex(ntuple(j -> j == i ? 1 : 0, N)), N)

    # 5. Main time-stepping loop
    for _ in 1:nsteps
        for I in inner_range
            laplacian = 0.0
            for i in 1:N
                # @inbounds removes safety checks for maximum speed
                @inbounds laplacian += r[i] * (u[I+e[i]] - 2.0 * u[I] + u[I-e[i]])
            end
            @inbounds u_new[I] = u[I] + laplacian
        end

        # O(1) pointer swap instead of copying data (u .= u_new)
        u, u_new = u_new, u
    end

    return PDESolution(axes_coords, u, p.T, p)
end


# ── Example Usage ─────────────────────────────────────────────────────────────

# 1D Problem
# prob1d = HeatEquation(
#     N_grid=(200,),
#     T=1.0,
#     f_init=x -> exp(-100 * (x - 0.5)^2)
# )
# sol1d = solve(prob1d)

# # 2D Problem
# prob2d = HeatEquation(
#     N_grid=(50, 50),
#     T=0.5,
#     f_init=(x, y) -> exp(-50 * ((x - 0.5)^2 + (y - 0.5)^2))
# )
# sol2d = solve(prob2d)

# # 3D Problem
# prob3d = HeatEquation(
#     N_grid=(20, 20, 20),
#     T=0.2,
#     f_init=(x, y, z) -> exp(-50 * ((x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2))
# )