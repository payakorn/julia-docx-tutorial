using Pkg
Pkg.activate(".")
using BenchmarkTools
using juliaPDEs

# ====================================================================
# 1. Define the "Slow" version (Without type parameters N and F)
# ====================================================================
Base.@kwdef struct HeatEquationSlow
    N_grid::Tuple              # Abstract! Julia doesn't know the size.
    Nt::Int    = 1000          # Number of time steps (Δt = T / Nt)
    α::Float64 = 0.01
    a::Tuple   = ntuple(_ -> 0.0, length(N_grid))   # Abstract Tuple
    b::Tuple   = ntuple(_ -> 1.0, length(N_grid))   # Abstract Tuple
    T::Float64 = 1.0
    f_init::Function           # Abstract! Julia doesn't know WHICH function.
    N::Int64   = length(N_grid)
end

function solve_slow(p::HeatEquationSlow)
    # N = length(p.N_grid) # We have to find N at runtime
    N = p.N

    d = ntuple(i -> (p.b[i] - p.a[i]) / (p.N_grid[i] + 1), N)
    axes_coords = ntuple(i -> collect(range(p.a[i] + d[i], p.b[i] - d[i], length=p.N_grid[i])), N)

    dt = p.T / p.Nt
    nsteps = p.Nt
    r = ntuple(i -> p.α * dt / d[i]^2, N)

    # Because N_grid is abstract, `zeros` creates a dynamically typed Array here
    u = zeros(Float64, p.N_grid)
    for I in CartesianIndices(u)
        coords = ntuple(dim -> axes_coords[dim][I[dim]], N)
        # 🚨 DYNAMIC DISPATCH:
        # Julia must pause and figure out what `f_init` actually is for every single iteration.
        u[I] = p.f_init(coords...)
    end

    for i in 1:N
        selectdim(u, i, 1) .= 0.0
        selectdim(u, i, p.N_grid[i]) .= 0.0
    end

    u_new = copy(u)
    inner_range = CartesianIndices(ntuple(i -> 2:(p.N_grid[i]-1), N))
    e = ntuple(i -> CartesianIndex(ntuple(j -> j == i ? 1 : 0, N)), N)

    # 🚨 TYPE INSTABILITY:
    # The lack of strict types cascades into the core simulation loop.
    for _ in 1:nsteps
        for I in inner_range
            laplacian = 0.0
            for i in 1:N
                laplacian += r[i] * (u[I+e[i]] - 2.0 * u[I] + u[I-e[i]])
            end
            u_new[I] = u[I] + laplacian
        end
        u, u_new = u_new, u
    end

    return u
end


# ====================================================================
# 2. Setup the test cases
# ====================================================================
println("Initializing 3D Heat Equation problems...")

N_pts = (70, 70, 70)
T_end = 0.5
Nt    = 5000          # number of time steps (same for both versions for a fair comparison)

# Fast version using the parameterized struct from juliaPDEs
fast_prob = HeatEquation(
    N_grid = N_pts,
    Nt     = Nt,
    T      = T_end,
    f_init = (x, y, z) -> exp(-100 * ((x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2))
)

# Slow version using our abstractly typed struct
slow_prob = HeatEquationSlow(
    N_grid = N_pts,
    Nt     = Nt,
    T      = T_end,
    f_init = (x, y, z) -> exp(-100 * ((x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2))
)

println("Warming up JIT compiler...")
solve(fast_prob)
solve_slow(slow_prob)


# ====================================================================
# 3. Run Benchmark
# ====================================================================
println("\n=== Benchmarking Fast (Parametric) Struct ===")
println("struct HeatEquation{N, F} ...")
@btime solve($fast_prob)
t_fast = @belapsed solve($fast_prob)

println("\n=== Benchmarking Slow (Abstract) Struct ===")
println("struct HeatEquationSlow ...")
@btime solve_slow($slow_prob)
t_slow = @belapsed solve_slow($slow_prob)

speedup = t_slow / t_fast

println("\n====================================================================")
println("🚀 RESULTS:")
println("The abstractly typed version is $(round(speedup, digits=1)) times SLOWER!")
println("====================================================================")
println("Notice the massive difference in execution time and memory allocations due to Type Boxing.")
