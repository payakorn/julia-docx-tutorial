using Pkg
Pkg.activate(".")
using LinearAlgebra, SparseArrays, BenchmarkTools, Printf, Statistics

# Same problem as the docs: 2-D heat on a 50x50 interior grid (~52x52 endpoint),
# Crank–Nicolson factorised once, then solved 500 times.
# Each phase is benchmarked with BenchmarkTools so we report median + min + std
# rather than a single noisy @elapsed sample.
const Nx = 50
const Ny = 50
const n  = Nx * Ny                      # total interior unknowns = 2500
const dx = 1.0 / (Nx + 1)
const α  = 0.05
const dt = 0.001
const θ  = 0.5
const NSTEPS = 500

println("Problem: 2-D heat,  $Nx × $Ny interior grid → n = $n unknowns")
println("Crank–Nicolson:     factor A once, solve $NSTEPS times")
println("Timer:              BenchmarkTools — median / min / std over many samples")
println("=" ^ 78)

# ────────────────────────────────────────────────────────────────────
# Build the 1-D tridiagonal Laplacian in BOTH representations
# ────────────────────────────────────────────────────────────────────
function tridi_dense(N, h)
    T = zeros(N, N)
    for i in 1:N
        T[i, i] = -2 / h^2
        i > 1 && (T[i, i-1] = 1 / h^2)
        i < N && (T[i, i+1] = 1 / h^2)
    end
    return T
end

tridi_sparse(N, h) = spdiagm(-1 => fill(1/h^2, N-1),
                              0 => fill(-2/h^2, N),
                              1 => fill(1/h^2, N-1))

# Small report helper — prints "median ± std (min)" in milliseconds.
function ms(trial::BenchmarkTools.Trial)
    med = median(trial).time / 1e6        # ns → ms
    mn  = minimum(trial).time / 1e6
    sd  = std(trial.times) / 1e6
    nsamples = length(trial.times)
    @sprintf("median %.3f ms  ±%.3f  (min %.3f, n=%d)", med, sd, mn, nsamples)
end

# Build the operators once for the matvec/solve benchmarks; rebuild them
# inside @benchmark for the assemble/factor benchmarks (each sample must
# do the work).
T1d_global = tridi_dense(Nx, dx)
T1s_global = tridi_sparse(Nx, dx)
Id_global  = Matrix(I, Nx, Nx)
Is_global  = sparse(I, Nx, Nx)
Ld_global  = kron(Id_global, T1d_global) .+ kron(T1d_global, Id_global)
Ls_global  = kron(Is_global, T1s_global) .+ kron(T1s_global, Is_global)
Ad_global  = Matrix(I, n, n)   .- (α * dt * θ) .* Ld_global
As_global  = sparse(I, n, n)   .- (α * dt * θ) .* Ls_global
Fd_global  = lu(Ad_global)
Fs_global  = lu(As_global)
b_global   = randn(n)

# ────────────────────────────────────────────────────────────────────
# 1. Assemble the 2-D Laplacian (rebuilt each sample)
# ────────────────────────────────────────────────────────────────────
println("\n[1] Assemble L = I⊗T + T⊗I  (Kron sum)")

t_assemble_dense = @benchmark begin
    T1 = tridi_dense($Nx, $dx)
    Id = Matrix(I, $Nx, $Nx)
    kron(Id, T1) .+ kron(T1, Id)
end seconds = 2
mem_dense_bytes = sizeof(Ld_global)
println("    dense    : ", ms(t_assemble_dense), "   |  memory ", mem_dense_bytes ÷ 1024^2, " MB")

t_assemble_sparse = @benchmark begin
    T1 = tridi_sparse($Nx, $dx)
    Is = sparse(I, $Nx, $Nx)
    kron(Is, T1) .+ kron(T1, Is)
end seconds = 2
nnz_sparse = nnz(Ls_global)
mem_sparse_bytes = Base.summarysize(Ls_global)
println("    sparse   : ", ms(t_assemble_sparse), "   |  memory ", mem_sparse_bytes ÷ 1024, " KB  |  nnz = ", nnz_sparse,
        @sprintf("  (%.3f%% of %d)", 100*nnz_sparse/n^2, n^2))

# ────────────────────────────────────────────────────────────────────
# 2. Build A = I − (α·dt·θ) L  and factor it
# ────────────────────────────────────────────────────────────────────
println("\n[2] One-time work — build A then lu(A)")

t_build_dense = @benchmark Matrix(I, $n, $n) .- ($α * $dt * $θ) .* $Ld_global  seconds = 1
t_build_sparse = @benchmark sparse(I, $n, $n) .- ($α * $dt * $θ) .* $Ls_global seconds = 1
println("    build A  dense : ", ms(t_build_dense))
println("    build A  sparse: ", ms(t_build_sparse))

t_factor_dense = @benchmark lu($Ad_global)  seconds = 2
t_factor_sparse = @benchmark lu($As_global) seconds = 2
println("    lu(A)    dense : ", ms(t_factor_dense))
println("    lu(A)    sparse: ", ms(t_factor_sparse))

# ────────────────────────────────────────────────────────────────────
# 3. CN time loop — per-step solve (Fd \ b) sampled many times
# ────────────────────────────────────────────────────────────────────
println("\n[3] One CN step:  F \\ b  (factor reused)")

t_solve_dense = @benchmark $Fd_global \ $b_global  seconds = 2
t_solve_sparse = @benchmark $Fs_global \ $b_global seconds = 2
println("    dense : ", ms(t_solve_dense))
println("    sparse: ", ms(t_solve_sparse))

# 500 steps end-to-end (single timed loop to capture loop overhead too)
println("\n[4] $NSTEPS Crank–Nicolson solves end-to-end")
t_loop_dense = @benchmark begin
    for _ in 1:$NSTEPS; $Fd_global \ $b_global; end
end seconds = 4
t_loop_sparse = @benchmark begin
    for _ in 1:$NSTEPS; $Fs_global \ $b_global; end
end seconds = 4
println("    dense : ", ms(t_loop_dense))
println("    sparse: ", ms(t_loop_sparse))

# ────────────────────────────────────────────────────────────────────
# Summary
# ────────────────────────────────────────────────────────────────────
println("\n" * "=" ^ 78)
println("SUMMARY (Nx=Ny=$Nx, n=$n unknowns, $NSTEPS CN solves)")
println("=" ^ 78)

mem_speedup    = mem_dense_bytes / mem_sparse_bytes
factor_speedup = median(t_factor_dense).time / median(t_factor_sparse).time
solve_speedup  = median(t_solve_dense).time / median(t_solve_sparse).time
loop_speedup   = median(t_loop_dense).time  / median(t_loop_sparse).time

med_ms(t) = median(t).time / 1e6
std_ms(t) = std(t.times)   / 1e6

@printf "  matrix memory   %8.1f MB  →  %8.1f KB         (×%5.0f smaller)\n" (mem_dense_bytes/1024^2) (mem_sparse_bytes/1024) mem_speedup
@printf "  factorisation   %7.2f ms ±%-5.2f  →  %7.3f ms ±%-5.3f  (×%4.0f faster)\n" med_ms(t_factor_dense) std_ms(t_factor_dense) med_ms(t_factor_sparse) std_ms(t_factor_sparse) factor_speedup
@printf "  one solve       %7.3f ms ±%-5.3f  →  %7.4f ms ±%-5.4f  (×%4.0f faster)\n" med_ms(t_solve_dense) std_ms(t_solve_dense) med_ms(t_solve_sparse) std_ms(t_solve_sparse) solve_speedup
@printf "  %d solves      %7.2f ms ±%-5.2f  →  %7.3f ms ±%-5.3f  (×%4.0f faster)\n" NSTEPS med_ms(t_loop_dense) std_ms(t_loop_dense) med_ms(t_loop_sparse) std_ms(t_loop_sparse) loop_speedup
