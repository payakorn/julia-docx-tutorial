using SparseArrays, LinearAlgebra

# ── 2D ────────────────────────────────────────────────────────────────────────
Base.@kwdef struct PoissonEquation <: EllipticProblem
    N::Int     = 50
    L::Float64 = 1.0
    # RHS forcing function f in  −∇²u = f(x,y)
    f::Function       = (x, y) -> 2π^2 * sin(π*x) * sin(π*y)
    # Known exact solution (nothing if not available)
    u_exact::Union{Function, Nothing} = (x, y) -> sin(π*x) * sin(π*y)
end

# ── Solver ─────────────────────────────────────────────────────────────────────
#
#  Discretises  −∇²u = f  on [0,L]² with Dirichlet BC u=0 on ∂Ω.
#
#  Grid:  N interior points per side, spacing h = L/(N+1).
#         x_i = i·h,  y_j = j·h,  i,j = 1…N
#
#  5-point stencil (scaled by h²):
#         −u_{i−1,j} − u_{i+1,j} − u_{i,j−1} − u_{i,j+1} + 4·u_{i,j} = h²·f_{ij}
#
#  In matrix form:  A · vec(U) = vec(F)
#  where  A = I_N ⊗ T₁ + T₁ ⊗ I_N   (Kronecker sum of 1-D Laplacians)
#  and vec() stacks columns (Julia column-major order).
#
function solve(p::PoissonEquation)
    N, h = p.N, p.L / (p.N + 1)

    x = collect(range(h, p.L - h, length=N))
    y = collect(range(h, p.L - h, length=N))

    # 1-D tridiagonal Laplacian  T₁ = h⁻² · tridiag(−1, 2, −1)
    T1 = spdiagm(
        -1 => fill(-1.0 / h^2, N-1),
         0 => fill( 2.0 / h^2, N),
         1 => fill(-1.0 / h^2, N-1)
    )
    IN = sparse(I, N, N)

    # 2-D Laplacian via Kronecker sum (column-major: x varies fastest)
    A = kron(IN, T1) + kron(T1, IN)      # N² × N² sparse, symmetric, positive definite

    # RHS: F[i,j] = f(x_i, y_j) — column-major so i is the fast index
    F   = [p.f(x[i], y[j]) for i in 1:N, j in 1:N]
    rhs = vec(F)

    # Direct sparse solve (uses UMFPACK under the hood)
    u_vec = A \ rhs
    u_mat = reshape(u_vec, N, N)

    return PDESolution((x, y), u_mat, 0.0, p)
end

# ── Error analysis helpers ────────────────────────────────────────────────────

"""
    l2_error(sol::PDESolution)

Compute the discrete L2 error against the exact solution stored in
`sol.problem.u_exact`. Returns `nothing` if no exact solution is available.
"""
function l2_error(sol::PDESolution{T, 2}) where T
    p = sol.problem
    p isa PoissonEquation || error("l2_error only defined for PoissonEquation")
    isnothing(p.u_exact) && return nothing

    N = size(sol, 1)
    h = p.L / (N + 1)

    err2 = 0.0
    for j in 1:N, i in 1:N
        e = sol[i, j] - p.u_exact(sol.x[i], sol.y[j])
        err2 += e^2
    end
    return h * sqrt(err2)          # discrete L2 norm ≈ ‖e‖_{L²}
end

"""
    convergence_table(Ns)

Run `PoissonEquation(N=n)` for each n in `Ns`, collect the L2 error, and
print a convergence table showing the observed order of accuracy.
"""
function convergence_table(Ns=[10, 20, 40, 80, 160])
    println(rpad("N", 6), rpad("h", 12), rpad("L² error", 14), "order")
    println("-"^40)
    prev_err = NaN
    for N in Ns
        sol = solve(PoissonEquation(N=N))
        h   = sol.problem.L / (N + 1)
        err = l2_error(sol)
        ord = isnan(prev_err) ? "—" : string(round(log2(prev_err / err), digits=2))
        println(rpad(N, 6),
                rpad(round(h,   digits=6),   12),
                rpad(round(err, sigdigits=4), 14),
                ord)
        prev_err = err
    end
end

# ── 3D ────────────────────────────────────────────────────────────────────────
#
#  −∇²u = f  on [0,L]³ with Dirichlet BCs u=0 on ∂Ω.
#
#  7-point stencil: 6·u_{i,j,k} − u_{i±1,j,k} − u_{i,j±1,k} − u_{i,j,k±1} = h²·f_{i,j,k}
#
#  Matrix: A = I⊗I⊗T1 + I⊗T1⊗I + T1⊗I⊗I  (3D Kronecker sum, N³×N³ sparse)
#
Base.@kwdef struct PoissonEquation3D <: EllipticProblem
    N::Int = 20
    L::Float64 = 1.0
    f::Function       = (x, y, z) -> 3π^2 * sin(π*x) * sin(π*y) * sin(π*z)
    u_exact::Union{Function, Nothing} = (x, y, z) -> sin(π*x) * sin(π*y) * sin(π*z)
end

function solve(p::PoissonEquation3D)
    N, h = p.N, p.L / (p.N + 1)

    x = collect(range(h, p.L - h, length=N))
    y = collect(range(h, p.L - h, length=N))
    z = collect(range(h, p.L - h, length=N))

    T1 = spdiagm(
        -1 => fill(-1.0 / h^2, N-1),
         0 => fill( 2.0 / h^2, N),
         1 => fill(-1.0 / h^2, N-1)
    )
    IN = sparse(I, N, N)

    A = kron(IN, kron(IN, T1)) + kron(IN, kron(T1, IN)) + kron(T1, kron(IN, IN))

    F   = [p.f(x[i], y[j], z[k]) for i in 1:N, j in 1:N, k in 1:N]
    u_vec = A \ vec(F)

    return PDESolution((x, y, z), reshape(u_vec, N, N, N), 0.0, p)
end

function l2_error(sol::PDESolution{T, 3}) where T
    p = sol.problem
    p isa PoissonEquation3D || error("l2_error only defined for PoissonEquation3D")
    isnothing(p.u_exact) && return nothing

    N = size(sol, 1)
    h = p.L / (N + 1)

    err2 = 0.0
    for k in 1:N, j in 1:N, i in 1:N
        e = sol[i, j, k] - p.u_exact(sol.x[i], sol.y[j], sol.z[k])
        err2 += e^2
    end
    return h^(3/2) * sqrt(err2)
end
