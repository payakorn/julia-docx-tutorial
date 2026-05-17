using SparseArrays, LinearAlgebra

# ── Poisson equation — dimension-free ─────────────────────────────────────────
#
#   −Δu = f(x...)          on the Cartesian box  ∏ᵢ [aᵢ, bᵢ]    in ℝᴺ
#   u = 0  on ∂Ω
#
# Parameters:
#   N   — spatial dimension (inferred from N_grid)
#   F   — type of the forcing function (specialised per closure)
#   EF  — type of the exact-solution function (or `Nothing` when unavailable —
#         that lets `l2_error` resolve the "no exact solution" branch statically)
#
Base.@kwdef struct PoissonEquation{N, F, EF} <: EllipticProblem
    N_grid::NTuple{N, Int} = (50, 50)                            # interior points per axis (defaults to 2D)
    a::NTuple{N, Float64} = ntuple(_ -> 0.0, length(N_grid))     # lower bound per axis
    b::NTuple{N, Float64} = ntuple(_ -> 1.0, length(N_grid))     # upper bound per axis
    f::F        = (x, y) -> 2π^2 * sin(π*x) * sin(π*y)           # default is 2D — override for other N
    u_exact::EF = (x, y) -> sin(π*x) * sin(π*y)                  # known exact solution, or `nothing`
end

function solve(p::PoissonEquation{N, F, EF}) where {N, F, EF}
    # 1. Build the per-axis interior grids via the shared helper.
    grids = ntuple(i -> interior_grid(p.a[i], p.b[i], p.N_grid[i]), N)
    h           = ntuple(i -> grids[i][1], N)
    axes_coords = ntuple(i -> grids[i][2], N)

    # 2. Assemble the N-D Laplacian as a Kronecker sum of axis-wise operators.
    A = nd_laplacian(h, p.N_grid)

    # 3. Build the right-hand side by sampling f on the Cartesian product grid.
    F_arr = [p.f(ntuple(dim -> axes_coords[dim][I[dim]], N)...) for I in CartesianIndices(p.N_grid)]

    # 4. Direct sparse solve (UMFPACK under the hood) and reshape back to N-D.
    u_vec = A \ vec(F_arr)
    u_arr = reshape(u_vec, p.N_grid)

    # 5. Wrap the solution with its coordinate grid.
    return PDESolution(axes_coords, u_arr, 0.0, p)
end

# ── N-D Laplacian via Kronecker sum ──────────────────────────────────────────
#
#  For each axis i, the discrete Laplacian acting along that axis is
#       Aᵢ = I_N ⊗ ⋯ ⊗ I_{i+1} ⊗ Tᵢ ⊗ I_{i-1} ⊗ ⋯ ⊗ I_1
#  where Tᵢ is the 1-D tridiagonal Laplacian h⁻² · tridiag(−1, 2, −1).
#  The full Laplacian is the sum A = Σᵢ Aᵢ.
#
#  Build the chain from the innermost (axis 1, fastest, rightmost in kron) to
#  the outermost (axis N, slowest, leftmost) — matches Julia's column-major
#  ordering convention.
#
function nd_laplacian(h::NTuple{N, Float64}, sizes::NTuple{N, Int}) where N
    total = prod(sizes)
    A     = spzeros(total, total)
    for axis in 1:N
        n  = sizes[axis]
        T1 = spdiagm(-1 => fill(-1.0 / h[axis]^2, n - 1),
                      0 => fill( 2.0 / h[axis]^2, n),
                      1 => fill(-1.0 / h[axis]^2, n - 1))
        op = T1
        for j in axis-1:-1:1              # wrap with identities on the inner (right) side
            op = kron(op, sparse(I, sizes[j], sizes[j]))
        end
        for j in axis+1:N                  # wrap with identities on the outer (left) side
            op = kron(sparse(I, sizes[j], sizes[j]), op)
        end
        A = A + op
    end
    return A
end

# ── Error analysis helpers ────────────────────────────────────────────────────

"""
    l2_error(sol::PDESolution)

Compute the discrete L² error against the exact solution stored in
`sol.problem.u_exact`. Returns `nothing` if no exact solution is available.

Works for any dimension N: the scaling is √(h₁·h₂·⋯·hₙ).
"""
function l2_error(sol::PDESolution{T, N}) where {T, N}
    p = sol.problem
    p isa PoissonEquation || error("l2_error only defined for PoissonEquation")
    isnothing(p.u_exact) && return nothing

    # 1. Recover the per-axis spacings from the same helper the solver used.
    h = ntuple(i -> interior_grid(p.a[i], p.b[i], p.N_grid[i])[1], N)

    # 2. Sum the squared point-wise errors over every grid node.
    err2 = 0.0
    for I in CartesianIndices(sol.u)
        coords = ntuple(dim -> sol.grid[dim][I[dim]], N)
        e = sol.u[I] - p.u_exact(coords...)
        err2 += e^2
    end

    # 3. Discrete L² norm scaling: √(h₁·…·hₙ).
    return sqrt(prod(h)) * sqrt(err2)
end

"""
    convergence_table(Ns; dims=2)

Run the canonical sin·sin… Poisson test problem at `dims` dimensions for each
`N` in `Ns`, collect the L² error, and print the empirical order column
`log₂(e_{N/2} / e_N)`. Should land on ~2.0 for the 5-point / 7-point stencil.
"""
function convergence_table(Ns=[10, 20, 40, 80, 160]; dims::Int=2)
    println(rpad("N", 6), rpad("h", 12), rpad("L² error", 14), "order")
    println("-"^40)
    prev_err = NaN
    for N in Ns
        prob = poisson_test_problem(N, dims)
        sol  = solve(prob)
        h, _ = interior_grid(prob.a[1], prob.b[1], N)
        err  = l2_error(sol)
        ord  = isnan(prev_err) ? "—" : string(round(log2(prev_err / err), digits=2))
        println(rpad(N, 6),
                rpad(round(h,   digits=6),   12),
                rpad(round(err, sigdigits=4), 14),
                ord)
        prev_err = err
    end
end

# Internal helper — the canonical sin-product test problem at any dimension.
function poisson_test_problem(N::Int, dims::Int)
    if dims == 2
        return PoissonEquation(
            N_grid  = (N, N),
            f       = (x, y) -> 2π^2 * sin(π*x) * sin(π*y),
            u_exact = (x, y) -> sin(π*x) * sin(π*y),
        )
    elseif dims == 3
        return PoissonEquation(
            N_grid  = (N, N, N),
            f       = (x, y, z) -> 3π^2 * sin(π*x) * sin(π*y) * sin(π*z),
            u_exact = (x, y, z) -> sin(π*x) * sin(π*y) * sin(π*z),
        )
    else
        error("poisson_test_problem only defined for dims = 2 or 3 (got $dims).")
    end
end
