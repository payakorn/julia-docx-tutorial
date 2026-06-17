# ---------------------------------------------------------------------------
# Makie tutorial examples
#
# A progression of small, self-contained figures showing how to build plots
# with CairoMakie — from a single line up to a multi-panel advanced figure.
# Every function returns a `Figure`; pass a path to `save(path, fig)` to write
# it to disk (e.g. `save("sin.png", makie_lines())`).
# ---------------------------------------------------------------------------

# 1. The simplest plot: one line --------------------------------------------
# `lines` is the workhorse for continuous curves. The first `Axis` argument
# can be a figure-grid slot (`f[1, 1]`) so Makie creates the axis for you.
function makie_simple_line()
  x = range(0, 2π, length=200)
  f = Figure()
  ax = Axis(f[1, 1], xlabel="x", ylabel="sin(x)", title="A single sine curve")
  CairoMakie.lines!(ax, x, sin.(x))
  return f
end

# 2. Several curves on one axis + a legend -----------------------------------
# Each `lines!` call adds another series; give it a `label` and call `axislegend`
# to draw the legend. `tan` is clipped with `ylims!` because it blows up.
function makie_sin_cos_tan()
  x = range(0, 2π, length=400)
  f = Figure()
  ax = Axis(f[1, 1], xlabel="x", ylabel="y", title="sin, cos and tan")

  CairoMakie.lines!(ax, x, sin.(x), label="sin", color=:steelblue)
  CairoMakie.lines!(ax, x, cos.(x), label="cos", color=:tomato)
  CairoMakie.lines!(ax, x, tan.(x), label="tan", color=:seagreen)

  CairoMakie.ylims!(ax, -3, 3)       # keep tan's asymptotes from dominating
  axislegend(ax; position=:rt)       # right-top corner
  return f
end

# 3. Markers and styling -----------------------------------------------------
# `scatter!` draws points; `lines!` with `linestyle` / `linewidth` styles a
# curve. Combine them to show data points over a fitted line.
function makie_scatter_line()
  x = range(0, 2π, length=25)
  noisy = sin.(x) .+ 0.15 .* randn(length(x))

  f = Figure()
  ax = Axis(f[1, 1], xlabel="x", ylabel="y", title="Data + reference curve")
  CairoMakie.lines!(ax, x, sin.(x), color=:gray, linestyle=:dash, label="sin(x)")
  CairoMakie.scatter!(ax, x, noisy, color=:purple, markersize=10, label="samples")
  axislegend(ax; position=:rb)
  return f
end

# 4. Subplots (a grid of axes) + a super-title -------------------------------
# `f[row, col]` places independent axes in a layout. `Label` spanning the top
# row (`f[0, :]`) acts as an overall title.
function makie_subplots()
  x = range(0, 2π, length=300)
  f = Figure(size=(800, 600))

  ax1 = Axis(f[1, 1], title="sin")
  ax2 = Axis(f[1, 2], title="cos")
  ax3 = Axis(f[2, 1], title="sin·cos")
  ax4 = Axis(f[2, 2], title="exp decay")

  CairoMakie.lines!(ax1, x, sin.(x), color=:steelblue)
  CairoMakie.lines!(ax2, x, cos.(x), color=:tomato)
  CairoMakie.lines!(ax3, x, sin.(x) .* cos.(x), color=:seagreen)
  CairoMakie.lines!(ax4, x, exp.(-x) .* sin.(4x), color=:darkorange)

  Label(f[0, :], "Four functions in a 2×2 grid", fontsize=20, font=:bold)
  return f
end

# 5. Heatmap + contour over a 2D field --------------------------------------
# Build a field z(x, y), then draw it as a filled `heatmap!` and overlay
# `contour!` lines. A shared `Colorbar` reads the heatmap's color scale.
function makie_contour()
  x = range(-3, 3, length=120)
  y = range(-3, 3, length=120)
  z = [exp(-(xi^2 + yi^2)) * sin(2xi) * cos(2yi) for xi in x, yi in y]

  f = Figure()
  ax = Axis(f[1, 1], xlabel="x", ylabel="y", title="Heatmap with contour lines",
    aspect=DataAspect())

  hm = CairoMakie.heatmap!(ax, x, y, z, colormap=:viridis)
  CairoMakie.contour!(ax, x, y, z, color=:white, linewidth=1, levels=10)
  Colorbar(f[1, 2], hm, label="z")
  return f
end

# 6. Filled contour (the original volcano example) ---------------------------
function makie_contourf()
  volcano = readdlm(Makie.assetpath("volcano.csv"), ',', Float64)

  f = Figure()
  ax = Axis(f[1, 1], title="Filled contour — Maunga Whau volcano")
  co = CairoMakie.contourf!(ax, volcano,
    levels=range(100, 180, length=10),
    extendlow=:cyan, extendhigh=:magenta)
  tightlimits!(ax)
  Colorbar(f[1, 2], co)
  return f
end

# 7. Advanced: 3D surface + a combined dashboard -----------------------------
# `Axis3` gives a perspective 3D axis for `surface!`. This figure also mixes
# axis types in one layout — a 3D surface beside a 2D heatmap — to show how
# different plot kinds compose in a single `Figure`.
function makie_advanced()
  x = range(-3, 3, length=80)
  y = range(-3, 3, length=80)
  z = [sin(xi) * cos(yi) * exp(-(xi^2 + yi^2) / 8) for xi in x, yi in y]

  f = Figure(size=(950, 450))

  ax1 = Axis3(f[1, 1], title="3D surface", xlabel="x", ylabel="y", zlabel="z",
    azimuth=0.6π, elevation=0.25π)
  CairoMakie.surface!(ax1, x, y, z, colormap=:plasma)

  ax2 = Axis(f[1, 2], title="Top-down heatmap", xlabel="x", ylabel="y",
    aspect=DataAspect())
  hm = CairoMakie.heatmap!(ax2, x, y, z, colormap=:plasma)
  CairoMakie.contour!(ax2, x, y, z, color=:black, linewidth=0.75, levels=8)
  Colorbar(f[1, 3], hm, label="z")

  Label(f[0, :], "Same field, two views", fontsize=20, font=:bold)
  return f
end

# ---------------------------------------------------------------------------
# Sharp output
#
# CairoMakie's PNG default is `px_per_unit=1`, which looks soft/blurry on hi-dpi
# screens and in print. Two cures:
#   * raster (PNG): raise `px_per_unit` (2–4) so more pixels back each figure
#     unit — crisp lines and text, larger file.
#   * vector (PDF/SVG/EPS): resolution-independent, stays razor-sharp at any
#     zoom; best for papers and slides.
#
# `save_sharp` picks sensible defaults: vector formats ignore `px_per_unit`
# (they don't need it), raster formats get a 3× scale unless you override it.
# ---------------------------------------------------------------------------
# `backend=CairoMakie` forces the file to be rendered through Cairo even when
# the GLMakie window backend is currently active — so saving a PDF/EPS keeps
# working no matter what you're using for on-screen display.
function save_sharp(path, fig; px_per_unit=3)
  ext = lowercase(splitext(path)[2])
  if ext in (".pdf", ".svg", ".eps")        # vector — already resolution-free
    Makie.save(path, fig; backend=CairoMakie)
  else                                       # raster — scale up the pixel grid
    Makie.save(path, fig; px_per_unit=px_per_unit, backend=CairoMakie)
  end
  return path
end

# Dump every tutorial figure to `dir` in the chosen format(s).
# e.g. save_makie_examples("figs", formats=("pdf", "png"))
function save_makie_examples(dir="makie_examples"; formats=("png",), px_per_unit=3)
  mkpath(dir)
  examples = [
    ("01_simple_line", makie_simple_line),
    ("02_sin_cos_tan", makie_sin_cos_tan),
    ("03_scatter_line", makie_scatter_line),
    ("04_subplots", makie_subplots),
    ("05_contour", makie_contour),
    ("06_contourf", makie_contourf),
    ("07_advanced", makie_advanced),
  ]
  paths = String[]
  for (name, fn) in examples
    fig = fn()
    for fmt in formats
      push!(paths, save_sharp(joinpath(dir, "$(name).$(fmt)"), fig; px_per_unit))
    end
  end
  return paths
end

# ---------------------------------------------------------------------------
# Choosing a display backend
#
# The example functions are backend-agnostic — `lines!`, `heatmap!`, `surface!`
# etc. are Makie-core functions, and whichever backend is *active* decides how a
# figure is shown. Two choices:
#
#   interactive()        -> GLMakie: a real GPU window, razor-sharp, zoom/rotate.
#                           Best for *looking at* plots (needs OpenGL; won't work
#                           over plain SSH / headless Docker).
#   sharp_display(n=3)   -> CairoMakie: static preview, rendered at n× pixel
#                           density so it doesn't look pixelated. Also the backend
#                           you want for saving PDF/EPS/SVG files.
#
# Switch any time; the same plotting code works under both.
# ---------------------------------------------------------------------------
function interactive()
  GLMakie.activate!()
  return :GLMakie
end

function sharp_display(px_per_unit::Real=3)
  CairoMakie.activate!(type="png", px_per_unit=px_per_unit)
  return px_per_unit
end

# Keep the original entry point working.
testmakie() = makie_contourf()
