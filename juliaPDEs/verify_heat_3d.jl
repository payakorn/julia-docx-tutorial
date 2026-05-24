using Pkg
Pkg.activate(".")
using juliaPDEs
using CairoMakie
using Printf, Statistics

# ── 3-D Heat verification ───────────────────────────────────────────
#
#     ∂u/∂t = α (u_xx + u_yy + u_zz),     (x, y, z) ∈ [0,1]³
#     u = 0 on every face,  u(x, y, z, 0) = sin(πx) sin(πy) sin(πz)
#
# The IC is a Laplacian eigenfunction (eigenvalue −3π²), so:
#
#     u(x, y, z, t) = sin(πx) sin(πy) sin(πz) · exp(−α · 3π² · t)
#
# Numerical solver: Crank–Nicolson (so temporal error doesn't pollute
# the spatial-convergence study). Memory grows as N³ so we stop the
# convergence study at N = 32 (≈ 33 000 unknowns, sparse LU fits in
# tens of MB and factors in a couple of seconds).
const α = 0.05
const T_END = 1.0
exact(x, y, z, t) = sin(π*x) * sin(π*y) * sin(π*z) * exp(-α * 3π^2 * t)

function solve_3d_cn(N::Int; T = T_END, Nt = 400)
    h    = 1.0 / N
    grid = Grid(a = (0.0, 0.0, 0.0), b = (1.0, 1.0, 1.0),
                stepsize = (h, h, h))
    zb   = (args...) -> 0.0
    tg   = TestGrid(grid,
                    (zb, zb, zb),
                    (zb, zb, zb),
                    (x, y, z) -> sin(π*x) * sin(π*y) * sin(π*z))
    return solve_implicit(HeatEquation(testgrid = tg, Nt = Nt, T = T, α = α);
                          θ = 0.5)
end

# ── Convergence study at fixed t = T_END ──────────────────────────
Ns      = [8, 12, 16, 24, 32]
errors  = Float64[]
hs      = Float64[]
println("3-D convergence at t = $T_END  (Crank–Nicolson, Nt = 400):")
for N in Ns
    t_start = time()
    sol = solve_3d_cn(N)
    ex  = [exact(x, y, z, T_END)
           for x in sol.x, y in sol.y, z in sol.z]
    e   = sol.u .- ex
    h   = sol.x[2] - sol.x[1]
    err = sqrt(sum(e .^ 2) * h^3)                # 3-D discrete L²
    push!(errors, err); push!(hs, h)
    @printf("  N = %3d   h = %.4e   ‖u_num − u_exact‖₂ = %.4e   (%.2fs)\n",
            N, h, err, time() - t_start)
end
println("\n  empirical orders:")
for k in 2:length(Ns)
    p = log(errors[k-1] / errors[k]) / log(hs[k-1] / hs[k])
    @printf("    N %3d → %3d   order ≈ %.3f\n", Ns[k-1], Ns[k], p)
end

# ── Mid-z slice for the visual panels ──────────────────────────────
const N_VIS = 40
sol_v   = solve_3d_cn(N_VIS)
mid     = div(size(sol_v.u, 3), 2)
z_mid   = sol_v.z[mid]
num_sl  = sol_v.u[:, :, mid]
ex_sl   = [exact(x, y, z_mid, T_END) for x in sol_v.x, y in sol_v.y]
err_sl  = num_sl .- ex_sl
err_mag = maximum(abs, err_sl)

println("\nVisualisation: $N_VIS × $N_VIS × $N_VIS, mid-plane z = $(round(z_mid, digits=3)),  max|error| = $(round(err_mag, sigdigits=3))")

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

fig = Figure(size = (1400, 1180))

# Shared scalar range for the heatmaps + the 3-D isosurfaces.
shared_lo = 0.0
shared_hi = max(maximum(sol_v.u), maximum([exact(x,y,z,T_END)
                 for x in sol_v.x, y in sol_v.y, z in sol_v.z]))

slice_label = @sprintf("z = %.3f", z_mid)

# Pre-build the exact field on the visualisation grid for the 3-D panel.
ex_full = [exact(x, y, z, T_END) for x in sol_v.x, y in sol_v.y, z in sol_v.z]

# 3-D viz: render a stack of z-slices as semi-transparent planes in 3-D
# space (CairoMakie has no real volume rendering but `surface!` with a
# constant-z height + a `color` matrix renders a flat coloured plane at
# that z — stack a handful of these and you get a "sandwich" view of
# the field that reveals the blob shape clearly).
#
# Plot back-to-front (highest z first) and use moderately high opacity
# so the bright middle layers aren't washed out by overlaying slices.
z_slice_vals = sort(collect(range(0.1, 0.9, length = 5)); rev = true)
function draw_slice_stack!(ax, field3d)
    for zv in z_slice_vals
        k     = argmin(abs.(sol_v.z .- zv))
        zfix  = sol_v.z[k]
        zsurf = fill(zfix, length(sol_v.x), length(sol_v.y))
        surface!(ax, sol_v.x, sol_v.y, zsurf;
                 color        = field3d[:, :, k],
                 colormap     = :thermal,
                 colorrange   = (shared_lo, shared_hi),
                 shading      = NoShading,
                 transparency = true,
                 alpha        = 0.82)
    end
end

# ── Row 1 — 3-D stacked-slice views ──────────────────────────────
ax3A = Axis3(fig[1, 1];
    title    = "Numerical  u_num(x, y, z, t = $T_END)  —  3-D stack",
    xlabel = "x", ylabel = "y", zlabel = "z",
    aspect = (1, 1, 1),
    azimuth = -π/4 + π/12,
    elevation = π/8,
    perspectiveness = 0.4,
    viewmode = :fitzoom,
)
draw_slice_stack!(ax3A, sol_v.u)

ax3B = Axis3(fig[1, 2];
    title    = "Exact  sin(πx) sin(πy) sin(πz) e^{−3απ²t}  —  3-D stack",
    xlabel = "x", ylabel = "y", zlabel = "z",
    aspect = (1, 1, 1),
    azimuth = -π/4 + π/12,
    elevation = π/8,
    perspectiveness = 0.4,
    viewmode = :fitzoom,
)
draw_slice_stack!(ax3B, ex_full)

# Shared colourbar for both 3-D panels
Colorbar(fig[1, 3];
    limits     = (shared_lo, shared_hi),
    colormap   = :thermal,
    label      = "u",
    width      = 14,
    height     = Relative(0.7),
)

# Caption row beneath the 3-D views
Label(fig[2, 1:3],
    "5 translucent z-slices at z = 0.1, 0.3, 0.5, 0.7, 0.9 (α = 0.82, back-to-front).  The bright middle slice (z = 0.5) carries the full sin(π z) peak; outer slices fade because sin(π · 0.1) = sin(π · 0.9) ≈ 0.31.";
    fontsize = 12, color = :gray30, halign = :center)
rowsize!(fig.layout, 2, Auto(0.04))

# ── Row 3 — mid-z slice heatmaps ──────────────────────────────────
axA = Axis(fig[3, 1], title = "Numerical slice  u_num(x, y, $slice_label)",
           xlabel = "x", ylabel = "y", aspect = DataAspect())
hmA = heatmap!(axA, sol_v.x, sol_v.y, num_sl,
               colormap = :thermal, colorrange = (shared_lo, shared_hi))

axB = Axis(fig[3, 2], title = "Exact slice  at  $slice_label",
           xlabel = "x", ylabel = "y", aspect = DataAspect())
hmB = heatmap!(axB, sol_v.x, sol_v.y, ex_sl,
               colormap = :thermal, colorrange = (shared_lo, shared_hi))
Colorbar(fig[3, 3], hmB, width = 12, height = Relative(0.85))

# ── Row 4 — error slice + convergence ─────────────────────────────
axC = Axis(fig[4, 1], title = "Pointwise error  u_num − u_exact  at  $slice_label",
           xlabel = "x", ylabel = "y", aspect = DataAspect())
hmC = heatmap!(axC, sol_v.x, sol_v.y, err_sl,
               colormap = :balance, colorrange = (-err_mag, err_mag))
Colorbar(fig[4, 2][1, 2], hmC, width = 12, height = Relative(0.85))

axD = Axis(fig[4, 2][1, 1],
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
text!(axD, Ns[1] * 1.3, errors[1] / 6;
      text = @sprintf("slope ≈ %.2f  ✓", order),
      align = (:left, :center),
      fontsize = 13, color = :crimson, font = :bold)
ylims!(axD, minimum(errors) / 4, maximum(errors) * 2)

# Row sizing — 3-D row taller than 2-D slice rows
rowsize!(fig.layout, 1, Relative(0.42))
rowsize!(fig.layout, 3, Relative(0.27))
rowsize!(fig.layout, 4, Relative(0.27))
colgap!(fig.layout, 24)
rowgap!(fig.layout, 18)

# ── Save ───────────────────────────────────────────────────────────
outdir  = joinpath(@__DIR__, "..", "site", "output", "figures")
mkpath(outdir)
outpath = joinpath(outdir, "heat_verification_3d.png")
save(outpath, fig; px_per_unit = 2)
println("\nSaved $outpath")
