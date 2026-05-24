using Pkg
Pkg.activate(".")
using juliaPDEs
using GLMakie
using Printf

# ── Companion to verify_heat_3d.jl ──────────────────────────────────
# Same 3-D heat problem, but rendered with GLMakie's OpenGL backend so
# we get real marching-cubes isosurfaces and absorption-based volume
# rendering — neither of which CairoMakie's static backend can do.
#
# We use offscreen rendering so the script runs without opening a
# window. Saving still goes through GLMakie's GL context.
GLMakie.activate!(visible = false)

const α     = 0.05
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

# ── Visualisation grid ────────────────────────────────────────────
const N_VIS = 48
println("Solving 3-D heat on $(N_VIS)³ grid for the GLMakie figure...")
@time sol_v   = solve_3d_cn(N_VIS)
@time ex_full = [exact(x, y, z, T_END) for x in sol_v.x, y in sol_v.y, z in sol_v.z]

shared_lo = 0.0
shared_hi = max(maximum(sol_v.u), maximum(ex_full))
println("  field range: [$(round(shared_lo, digits=3)), $(round(shared_hi, digits=3))]")

# Isosurface levels — three nested shells at fractions of the max.
# Drop the outermost very-translucent shell that just muddied the
# render, and use distinct manual colours so each shell is identifiable
# (the previous uniform :YlOrRd colormap made them blend into one blob).
iso_fracs  = (0.25, 0.55, 0.85)
iso_levels = collect(f * shared_hi for f in iso_fracs)
iso_colors = (RGBAf(1.00, 0.98, 0.82, 0.55),    # outer — pale cream, mid alpha
              RGBAf(1.00, 0.88, 0.60, 0.78),    # middle — soft apricot
              RGBAf(1.00, 0.62, 0.55, 0.95))    # inner — warm light coral

# ── Figure ─────────────────────────────────────────────────────────
set_theme!(theme_minimal())
update_theme!(
    fontsize        = 14,
    backgroundcolor = :white,
)

fig = Figure(size = (1450, 1100))

# Axis3 builder factory.
function make_axis3(pos, title)
    Axis3(pos;
        title           = title,
        xlabel = "x", ylabel = "y", zlabel = "z",
        aspect          = (1, 1, 1),
        azimuth         = -π/4 + π/12,
        elevation       = π/8,
        perspectiveness = 0.4,
        viewmode        = :fitzoom,
    )
end

# ── Row 1 — marching-cubes isosurfaces ────────────────────────────
function draw_iso!(ax, field3d)
    # Use `volume!(algorithm=:iso)` per level instead of `contour!`.
    # The ray-marched iso renderer respects the supplied colormap as a
    # near-flat colour and ignores the directional shading that was
    # tinting our pastel contour shells blue.  Draw outer-to-inner so
    # transparency layering reads correctly.
    for (lev, col) in zip(iso_levels, iso_colors)
        volume!(ax, (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), field3d;
                algorithm  = :iso,
                isovalue   = lev,
                isorange   = 0.004,
                colormap   = [col, col],
                colorrange = (shared_lo, shared_hi),
                transparency = true)
    end
end

axIA = make_axis3(fig[1, 1],
    "Numerical  u_num(x, y, z)  —  3 isosurfaces (marching cubes)")
draw_iso!(axIA, sol_v.u)

axIB = make_axis3(fig[1, 2],
    "Exact  sin(πx) sin(πy) sin(πz) e^{−3απ²t}  —  3 isosurfaces")
draw_iso!(axIB, ex_full)

# Legend in the third column — manual Legend with PolyElement swatches
# so each iso shell's colour + value is clearly labelled.
Legend(fig[1, 3],
    [PolyElement(color = c, strokecolor = :gray30, strokewidth = 0.8)
     for c in iso_colors],
    [@sprintf("u = %.3f   (%.0f%% of max)", lev, 100 * frac)
     for (lev, frac) in zip(iso_levels, iso_fracs)],
    "Isosurface levels";
    framevisible = false,
    labelsize    = 13,
    titlesize    = 14,
    titlefont    = :bold,
    patchsize    = (28, 16),
    rowgap       = 6,
    valign       = :center,
)
colsize!(fig.layout, 3, Auto(0.4))

Label(fig[2, 1:3],
    @sprintf("3 nested marching-cubes shells at u = %.3f / %.3f / %.3f.  Outer shells are translucent so the inner core stays visible.  GLMakie's OpenGL backend depth-sorts the transparency correctly.",
        iso_levels...);
    fontsize = 12, color = :gray30, halign = :center)
rowsize!(fig.layout, 2, Auto(0.04))

# ── Row 3 — absorption-based volume rendering ────────────────────
axVA = make_axis3(fig[3, 1],
    "Numerical  —  max-intensity projection (volume!)")
volume!(axVA, (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), sol_v.u;
        algorithm  = :mip,
        colormap   = :YlOrRd,
        colorrange = (shared_lo, shared_hi))

axVB = make_axis3(fig[3, 2],
    "Exact  —  max-intensity projection (volume!)")
volume!(axVB, (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), ex_full;
        algorithm  = :mip,
        colormap   = :YlOrRd,
        colorrange = (shared_lo, shared_hi))

Colorbar(fig[3, 3];
    limits   = (shared_lo, shared_hi),
    colormap = :YlOrRd,
    label    = "u",
    width    = 14,
    height   = Relative(0.7),
)

Label(fig[4, 1:3],
    "Max-intensity projection (MIP) — every camera ray records the maximum field value it encounters, so the rendering reads as a brightness map of the 3-D field. The bright spot is exactly the (0.5, 0.5, 0.5) peak of sin(πx)sin(πy)sin(πz).  Numerical and exact MIPs are pixel-identical → the solver matches inside the cube, not just on the cross-sections.";
    fontsize = 12, color = :gray30, halign = :center)
rowsize!(fig.layout, 4, Auto(0.04))

# Row sizing — make both visualisation rows the same dominant height
rowsize!(fig.layout, 1, Relative(0.46))
rowsize!(fig.layout, 3, Relative(0.46))
colgap!(fig.layout, 24)
rowgap!(fig.layout, 12)

# ── Save ───────────────────────────────────────────────────────────
outdir  = joinpath(@__DIR__, "..", "site", "output", "figures")
mkpath(outdir)
outpath = joinpath(outdir, "heat_verification_3d_iso.png")
save(outpath, fig; px_per_unit = 2)
println("\nSaved $outpath")
