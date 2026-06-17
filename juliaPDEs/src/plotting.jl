# ── Plotting a PDESolution with Makie ─────────────────────────────────────────
#
# `plot_solution` turns any `PDESolution` into a Makie `Figure`, dispatching on
# the spatial dimension N:
#
#   N = 1  →  line plot         u(x)
#   N = 2  →  heatmap (default), or :surface / :contour / :contourf
#   N = 3  →  stacked colored z-slices inside an Axis3
#
# Makie is used qualified (`Makie.*`) so its exports never clash with the `Plots`
# functions the rest of the package uses. To actually render or save the figure,
# load a backend in your session:
#
#   using CairoMakie        # static PNG/SVG, headless
#   fig = plot_solution(sol)
#   save("solution.png", fig)
#
import Makie

"""
    plot_solution(sol::PDESolution; kind = :auto, colormap = :viridis,
                  title = nothing, slices = (0.25, 0.5, 0.75)) -> Makie.Figure

Build a Makie figure for `sol`, choosing a sensible plot for its dimension.

- `kind`     — `:auto` (default), or for 2-D fields `:heatmap`, `:surface`,
               `:contour`, `:contourf`.
- `colormap` — any Makie colormap symbol.
- `slices`   — for 3-D fields, the fractional z-heights to draw as planes.
"""
function plot_solution(sol::PDESolution{T,N}; kind::Symbol = :auto,
                       colormap = :viridis, title = nothing,
                       slices = (0.25, 0.5, 0.75)) where {T,N}
    if N == 1
        return _plot_1d(sol; colormap, title)
    elseif N == 2
        return _plot_2d(sol; kind, colormap, title)
    elseif N == 3
        return _plot_3d_slices(sol; colormap, title, slices)
    else
        error("plot_solution: no plot defined for $(N)-D solutions.")
    end
end

_default_title(sol) = string(nameof(typeof(getfield(sol, :problem))),
                             "  (t = ", round(sol.t; digits=4), ")")

function _plot_1d(sol; colormap, title)
    fig = Makie.Figure(size = (760, 460), fontsize = 16)
    ax  = Makie.Axis(fig[1, 1]; title = something(title, _default_title(sol)),
                     xlabel = "x", ylabel = "u")
    Makie.lines!(ax, sol.x, sol.u; linewidth = 2.5)
    return fig
end

function _plot_2d(sol; kind, colormap, title)
    kind === :auto && (kind = :heatmap)
    ttl = something(title, _default_title(sol))

    if kind === :surface
        fig = Makie.Figure(size = (820, 640), fontsize = 16)
        ax  = Makie.Axis3(fig[1, 1]; title = ttl, xlabel = "x", ylabel = "y",
                          zlabel = "u", azimuth = 0.6π, elevation = 0.18π)
        sp  = Makie.surface!(ax, sol.x, sol.y, sol.u; colormap)
        Makie.Colorbar(fig[1, 2], sp; label = "u(x,y)")
        return fig
    end

    fig = Makie.Figure(size = (760, 620), fontsize = 16)
    ax  = Makie.Axis(fig[1, 1]; title = ttl, xlabel = "x", ylabel = "y",
                     aspect = Makie.DataAspect())
    if kind === :heatmap
        pl = Makie.heatmap!(ax, sol.x, sol.y, sol.u; colormap)
    elseif kind === :contourf || kind === :contour
        levels = range(minimum(sol.u), maximum(sol.u); length = 12)
        pl = Makie.contourf!(ax, sol.x, sol.y, sol.u; levels, colormap)
        kind === :contour && Makie.contour!(ax, sol.x, sol.y, sol.u;
                                            levels, color = :black, linewidth = 0.9)
    else
        error("plot_solution: unknown kind=$(repr(kind)) for a 2-D solution " *
              "(use :heatmap, :surface, :contour, or :contourf).")
    end
    Makie.Colorbar(fig[1, 2], pl; label = "u(x,y)")
    return fig
end

function _plot_3d_slices(sol; colormap, title, slices)
    Nz  = length(sol.z)
    lo, hi = extrema(sol.u)
    fig = Makie.Figure(size = (820, 660), fontsize = 16)
    ax  = Makie.Axis3(fig[1, 1]; title = something(title, _default_title(sol)),
                      xlabel = "x", ylabel = "y", zlabel = "z",
                      azimuth = 0.55π, elevation = 0.16π)
    local pl
    for frac in slices
        k      = clamp(round(Int, frac * Nz), 1, Nz)
        zplane = fill(sol.z[k], length(sol.x), length(sol.y))
        pl = Makie.surface!(ax, sol.x, sol.y, zplane;
                            color = sol.u[:, :, k], colormap,
                            colorrange = (lo, hi), transparency = true)
    end
    Makie.Colorbar(fig[1, 2], pl; label = "u(x,y,z)")
    return fig
end
