using Pkg
Pkg.activate(".")
using juliaPDEs
using CairoMakie
using Printf, Statistics

# ── 2-D Heat verification ───────────────────────────────────────────
#
#     ∂u/∂t = α (∂²u/∂x² + ∂²u/∂y²),     (x, y) ∈ [0,1]²
#     u = 0 on every face,  u(x, y, 0) = sin(πx) sin(πy)
#
# The IC is a Laplacian eigenfunction (eigenvalue −2π²), so the exact
# solution is one decaying mode:
#
#     u(x, y, t) = sin(πx) sin(πy) · exp(−α · 2π² · t)
#
# Numerical solver: Crank–Nicolson (so temporal error doesn't pollute
# the spatial-convergence study at small Δx).
const α = 0.05
const T_END = 1.0
exact(x, y, t) = sin(π * x) * sin(π * y) * exp(-α * 2π^2 * t)

function solve_2d_cn(N::Int; T = T_END, Nt = 400)
    h    = 1.0 / N
    grid = Grid(a = (0.0, 0.0), b = (1.0, 1.0), stepsize = (h, h))
    zb   = (args...) -> 0.0
    tg   = TestGrid(grid,
                    (zb, zb),
                    (zb, zb),
                    (x, y) -> sin(π * x) * sin(π * y))
    return solve_implicit(HeatEquation(testgrid = tg, Nt = Nt, T = T, α = α);
                          θ = 0.5)
end

# ── Convergence study at fixed t = T_END ──────────────────────────
Ns      = [16, 32, 64, 128]
errors  = Float64[]
hs      = Float64[]
println("2-D convergence at t = $T_END  (Crank–Nicolson, Nt = 400):")
for N in Ns
    sol = solve_2d_cn(N)
    ex  = [exact(x, y, T_END) for x in sol.x, y in sol.y]
    e   = sol.u .- ex
    h   = sol.x[2] - sol.x[1]
    err = sqrt(sum(e .^ 2) * h^2)                # 2-D discrete L²
    push!(errors, err); push!(hs, h)
    @printf("  N = %4d   h = %.4e   ‖u_num − u_exact‖₂ = %.4e\n", N, h, err)
end
println("\n  empirical orders:")
for k in 2:length(Ns)
    p = log(errors[k-1] / errors[k]) / log(hs[k-1] / hs[k])
    @printf("    N %4d → %4d   order ≈ %.3f\n", Ns[k-1], Ns[k], p)
end

# ── Heatmap-resolution solve for the visual panels ────────────────
const N_VIS = 80
sol_v = solve_2d_cn(N_VIS)
ex_v  = [exact(x, y, T_END) for x in sol_v.x, y in sol_v.y]
err_v = sol_v.u .- ex_v
err_mag = maximum(abs, err_v)

println("\nVisualisation grid: $N_VIS × $N_VIS,  max|error| = $(round(err_mag, sigdigits=3))")

# ── Figure ─────────────────────────────────────────────────────────
set_theme!(theme_minimal())
update_theme!(
    fontsize        = 14,
    backgroundcolor = :white,
    Axis = (
        backgroundcolor   = :white,
        xgridcolor        = (:gray, 0.18),
        ygridcolor        = (:gray, 0.18),
        spinewidth        = 1.0,
        xtickwidth        = 1.0, ytickwidth = 1.0,
        titlefont         = :regular,
        titlesize         = 15,
        titlegap          = 8,
        xlabelfont        = :regular, ylabelfont = :regular,
        xlabelsize        = 13, ylabelsize = 13,
    ),
)

fig = Figure(size = (1280, 800))

# Shared colour range for numerical & exact so the visual comparison is fair.
shared_lo = 0.0
shared_hi = max(maximum(sol_v.u), maximum(ex_v))

# Panel A — numerical
axA = Axis(fig[1, 1], title = "Numerical  u_num(x, y, t = $T_END)",
           xlabel = "x", ylabel = "y", aspect = DataAspect())
hmA = heatmap!(axA, sol_v.x, sol_v.y, sol_v.u,
               colormap = :thermal, colorrange = (shared_lo, shared_hi))
Colorbar(fig[1, 2], hmA, width = 12)

# Panel B — exact
axB = Axis(fig[1, 3], title = "Exact  sin(πx) sin(πy) e^{−2απ²t}",
           xlabel = "x", ylabel = "y", aspect = DataAspect())
hmB = heatmap!(axB, sol_v.x, sol_v.y, ex_v,
               colormap = :thermal, colorrange = (shared_lo, shared_hi))
Colorbar(fig[1, 4], hmB, width = 12)

# Panel C — pointwise error (diverging map, symmetric range)
axC = Axis(fig[2, 1], title = "Pointwise error  u_num − u_exact",
           xlabel = "x", ylabel = "y", aspect = DataAspect())
hmC = heatmap!(axC, sol_v.x, sol_v.y, err_v,
               colormap = :balance, colorrange = (-err_mag, err_mag))
Colorbar(fig[2, 2], hmC, width = 12)

# Panel D — convergence
axD = Axis(fig[2, 3:4],
           title  = "Convergence — second-order in Δx",
           xlabel = "grid points per axis  N",
           ylabel = "L² error at t = $T_END",
           xscale = log10, yscale = log10,
           xticks = (Ns, string.(Ns)))
scatterlines!(axD, Ns, errors,
              color = :crimson, marker = :diamond, markersize = 12,
              linewidth = 2.2, label = "measured")
ref = errors[1] .* (hs ./ hs[1]) .^ 2
lines!(axD, Ns, ref, color = :gray40, linestyle = :dash, linewidth = 1.6,
       label = "O(Δx²) reference")
axislegend(axD, position = :lb, framevisible = false, labelsize = 12)

order = log(errors[end] / errors[1]) / log(hs[end] / hs[1])
text!(axD, 25, errors[1] / 6;
      text = @sprintf("slope ≈ %.2f  ✓", order),
      align = (:left, :center),
      fontsize = 13, color = :crimson, font = :bold)
ylims!(axD, minimum(errors) / 4, maximum(errors) * 2)

# ── Save ───────────────────────────────────────────────────────────
outdir  = joinpath(@__DIR__, "..", "site", "output", "figures")
mkpath(outdir)
outpath = joinpath(outdir, "heat_verification_2d.png")
save(outpath, fig; px_per_unit = 2)
println("\nSaved $outpath")
