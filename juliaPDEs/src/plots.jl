# ── Figure functions — dispatch on Problem structs ────────────────────────────
#
# Each `fig_*` function accepts a concrete Problem struct, runs the solver
# (internally capturing snapshot history where needed), and returns a
# `Plots.Plot` object. Pass `savepath` to also write to disk.
#
# `plot_solution(sol::PDESolution; kwargs...)` dispatches on `sol.problem`.
#
using Printf

# ══════════════════════════════════════════════════════════════════════════════
# Internal: 1-D heat — forward Euler with snapshot + space-time history
# ══════════════════════════════════════════════════════════════════════════════

function _heat_1d_snapshots(p::HeatEquation{1,F};
        snapshot_times::Vector{Float64} = [0.0, 0.05, 0.2, 0.5, 1.0],
        save_every::Int = 0) where F
    tg    = p.testgrid
    x     = tg.grid.space[1]
    npts  = tg.grid.numgrid[1]
    dt    = p.T / p.Nt
    dx    = (x[end] - x[1]) / (npts - 1)
    r     = p.α * dt / dx^2

    u      = [tg.f_init(x[i]) for i in 1:npts]
    u[1]   = tg.bc1[1](0.0)
    u[end] = tg.bc2[1](0.0)

    snaps      = Dict{Float64, Vector{Float64}}()
    snaps[0.0] = copy(u)
    targets    = sort(filter(t -> t > 0.0, snapshot_times))

    n_save = save_every > 0 ? save_every : max(1, p.Nt ÷ 200)
    hist   = [copy(u)]
    t_hist = [0.0]

    for step in 1:p.Nt
        t = step * dt
        u[2:end-1] .+= r .* (u[3:end] .- 2 .* u[2:end-1] .+ u[1:end-2])
        u[1]   = tg.bc1[1](t)
        u[end] = tg.bc2[1](t)
        for tt in targets
            !haskey(snaps, tt) && t >= tt && (snaps[tt] = copy(u))
        end
        step % n_save == 0 && (push!(hist, copy(u)); push!(t_hist, t))
    end
    return x, snaps, hcat(hist...), t_hist
end

# ══════════════════════════════════════════════════════════════════════════════
# Internal: 1-D wave — leap-frog with snapshot + space-time history
# ══════════════════════════════════════════════════════════════════════════════

function _wave_1d_snapshots(p::WaveEquation{1,F};
        snapshot_times::Vector{Float64} = Float64[],
        save_every::Int = 0) where F
    dx, x = endpoint_grid(p.a[1], p.b[1], p.N_grid[1])
    N      = p.N_grid[1]
    dt     = p.T / p.Nt
    λ      = (p.c * dt / dx)^2

    u_prev    = [p.f_init(x[i]) for i in 1:N]
    u_prev[1] = u_prev[end] = 0.0
    u_curr    = copy(u_prev)
    u_next    = similar(u_curr)

    snaps      = Dict{Float64, Vector{Float64}}()
    snaps[0.0] = copy(u_curr)
    targets    = sort(filter(t -> t > 0.0, snapshot_times))

    n_save = save_every > 0 ? save_every : max(1, p.Nt ÷ 250)
    hist   = [copy(u_curr)]
    t_hist = [0.0]
    t      = 0.0

    for step in 1:p.Nt
        u_next[2:end-1] .= 2 .* u_curr[2:end-1] .- u_prev[2:end-1] .+
            λ .* (u_curr[3:end] .- 2 .* u_curr[2:end-1] .+ u_curr[1:end-2])
        u_next[1] = u_next[end] = 0.0
        u_prev, u_curr, u_next = u_curr, u_next, u_prev
        t += dt
        for tt in targets
            !haskey(snaps, tt) && t >= tt && (snaps[tt] = copy(u_curr))
        end
        step % n_save == 0 && (push!(hist, copy(u_curr)); push!(t_hist, t))
    end
    return x, snaps, hcat(hist...), t_hist
end

# ══════════════════════════════════════════════════════════════════════════════
# Public figure functions
# ══════════════════════════════════════════════════════════════════════════════

"""
    fig_heat_equation(prob::HeatEquation{1}; snapshot_times, savepath, size, dpi) -> Plot

2-panel figure: temperature snapshots at `snapshot_times` (left) and a full
space-time heatmap (right). Returns a `Plots.Plot`; writes to `savepath` if
non-empty.
"""
function fig_heat_equation(p::HeatEquation{1,F};
        snapshot_times::Vector{Float64} = [0.0, 0.05, 0.2, 0.5, 1.0],
        savepath::String = "",
        size::Tuple{Int,Int} = (1100, 450),
        dpi::Int = 150) where F

    x, snaps, U_hist, t_hist = _heat_1d_snapshots(p; snapshot_times)
    ts    = sort(collect(keys(snaps)))
    n     = length(ts)
    cmap  = cgrad(:YlOrRd, range(0.3, 0.95, length=n))
    ylo   = minimum(minimum(v) for v in values(snaps))
    yhi   = maximum(maximum(v) for v in values(snaps))
    pad   = max(0.05 * abs(yhi - ylo), 1e-3)

    p1 = plot(xlabel="x", ylabel="u(x,t)",
              title="Heat Equation:  ∂u/∂t = α∇²u   (α = $(p.α))",
              legend=:topright, grid=true, gridalpha=0.3,
              ylim=(ylo - pad, yhi + pad), framestyle=:box)
    for (i, tt) in enumerate(ts)
        plot!(p1, x, snaps[tt], lw=2.2, color=cmap[i],
              label=@sprintf("t = %.2f", tt))
    end
    plot!(p1, [x[end÷2], x[end÷2]], [yhi * 0.9, yhi * 0.25],
          arrow=true, color=:firebrick, lw=1.5, label="")
    annotate!(p1, x[end÷2] + 0.08, (yhi * 0.9 + yhi * 0.25) / 2,
              text("diffuses\n& smooths", :firebrick, :left, 9))

    p2 = heatmap(t_hist, x, U_hist, c=:YlOrRd,
                 xlabel="time t", ylabel="x",
                 title="Space-time evolution",
                 colorbar_title="u")

    plt = plot(p1, p2, layout=(1,2), size=size, dpi=dpi,
               plot_title="Heat equation — α = $(p.α),  T = $(p.T)",
               plot_titlefontsize=13)
    isempty(savepath) || savefig(plt, savepath)
    return plt
end

"""
    fig_wave_equation(prob::WaveEquation{1}; snapshot_times, T_spacetime, savepath, size, dpi) -> Plot

2-panel figure: displacement snapshots (left) and a space-time heatmap showing
characteristic lines (right). `T_spacetime` sets the final time for the longer
space-time panel (default: 2 × prob.T).
"""
function fig_wave_equation(p::WaveEquation{1,F};
        snapshot_times::Vector{Float64} = Float64[],
        T_spacetime::Float64 = 2.0 * p.T,
        savepath::String = "",
        size::Tuple{Int,Int} = (1100, 450),
        dpi::Int = 150) where F

    if isempty(snapshot_times)
        snapshot_times = sort(unique([0.0; [k * p.T / 4 for k in 1:4]]))
    end

    x, snaps, _, _ = _wave_1d_snapshots(p; snapshot_times)
    ts      = sort(collect(keys(snaps)))
    n       = length(ts)
    cmap    = cgrad(:viridis, range(0.1, 0.9, length=n))
    max_amp = maximum(abs.(p.f_init.(x)))        # IC amplitude ≈ conserved
    ylim_v  = (-(max_amp + 0.1), max_amp + 0.1)

    p1 = plot(xlabel="x", ylabel="u(x,t)",
              title="Wave Equation:  ∂²u/∂t² = c²∇²u   (c = $(p.c))",
              legend=:topright, grid=true, gridalpha=0.3,
              ylim=ylim_v, framestyle=:box)
    hline!(p1, [0], color=:gray, lw=0.5, label="")
    for (i, tt) in enumerate(ts)
        plot!(p1, x, snaps[tt], lw=2.0, color=cmap[i],
              label=@sprintf("t = %.2f", tt))
    end

    # Longer simulation for the space-time panel
    Nt_long = ceil(Int, p.Nt * T_spacetime / p.T)
    p_long  = WaveEquation(N_grid=p.N_grid, a=p.a, b=p.b,
                           Nt=Nt_long, T=T_spacetime,
                           c=p.c, f_init=p.f_init)
    _, _, U_long, t_long = _wave_1d_snapshots(p_long)

    p2 = heatmap(t_long, x, U_long, c=:RdBu, clims=(-max_amp, max_amp),
                 xlabel="time t", ylabel="x",
                 title="Space-time: characteristic lines",
                 colorbar_title="u")

    plt = plot(p1, p2, layout=(1,2), size=size, dpi=dpi,
               plot_title="Wave equation — c = $(p.c),  T = $(p.T)",
               plot_titlefontsize=13)
    isempty(savepath) || savefig(plt, savepath)
    return plt
end

"""
    fig_poisson_equation(prob::PoissonEquation{2}; savepath, size, dpi) -> Plot

3-panel figure: source term f(x,y), solution u(x,y), and a 3-D surface view.
"""
function fig_poisson_equation(p::PoissonEquation{2,F,EF};
        savepath::String = "",
        size::Tuple{Int,Int} = (1300, 480),
        dpi::Int = 150) where {F, EF}

    sol = solve(p)
    xs, ys = sol.x, sol.y
    # sol.u has shape (Nx, Ny): sol.u[i,j] = u(xs[i], ys[j])
    # heatmap(x, y, Z) expects Z[j,i] = Z(xs[i], ys[j]) → transpose
    Fm = [p.f(xs[i], ys[j]) for j in eachindex(ys), i in eachindex(xs)]

    p1 = heatmap(xs, ys, Fm, c=:coolwarm,
                 xlabel="x", ylabel="y",
                 title="Source term  f(x,y)", aspect_ratio=:equal)
    p2 = heatmap(xs, ys, sol.u', c=:viridis,
                 xlabel="x", ylabel="y",
                 title="Solution  u(x,y)", aspect_ratio=:equal)
    p3 = surface(xs, ys, sol.u', c=:viridis, camera=(40, 30),
                 xlabel="x", ylabel="y", zlabel="u",
                 title="3D view of u(x,y)")

    plt = plot(p1, p2, p3, layout=(1,3), size=size, dpi=dpi,
               plot_title="Poisson Equation:  −∇²u = f   (steady-state, no time)",
               plot_titlefontsize=14)
    isempty(savepath) || savefig(plt, savepath)
    return plt
end

"""
    fig_navier_stokes(prob::LidCavityFlow; savepath, size, dpi) -> Plot

2-panel figure: lid-driven cavity schematic (left) and a speed-magnitude
heatmap with ψ-streamlines overlaid (right).
"""
function fig_navier_stokes(p::LidCavityFlow;
        savepath::String = "",
        size::Tuple{Int,Int} = (1100, 500),
        dpi::Int = 150)

    ψ, _, _, _, speed, _ = _solve_psi_omega(p)
    N  = p.N
    xg = range(0, 1, length=N)
    yg = range(0, 1, length=N)
    θ  = range(0, 2π, length=100)

    p1 = plot(legend=false, framestyle=:none, aspect_ratio=:equal,
              xlim=(-0.2, 1.2), ylim=(-0.2, 1.4),
              title="Lid-Driven Cavity (NS benchmark)")
    plot!(p1, [0,1,1,0,0], [0,0,1,1,0], color=:black, lw=2.5)
    plot!(p1, [0.15, 0.85], [1.05, 1.05], arrow=true, color=:firebrick, lw=2.5)
    annotate!(p1,  0.5,  1.18, text("moving lid  U = $(p.U_lid)", :firebrick, :center, 11, :bold))
    annotate!(p1, -0.07, 0.5,  text("u=v=0", :steelblue, :right,   10))
    annotate!(p1,  1.07, 0.5,  text("u=v=0", :steelblue, :left,    10))
    annotate!(p1,  0.5, -0.10, text("u=v=0", :steelblue, :center,  10))
    for r in (0.15, 0.25, 0.35)
        plot!(p1, 0.5 .+ r .* cos.(θ), 0.55 .+ r .* sin.(θ),
              color=:seagreen, lw=1, alpha=0.6)
    end
    annotate!(p1, 0.5, 0.55, text("vortex", :seagreen, :center, 11, :bold))

    p2 = heatmap(xg, yg, speed', c=:plasma, aspect_ratio=:equal,
                 xlims=(0,1), ylims=(0,1),
                 xlabel="x", ylabel="y",
                 title=@sprintf("Numerical |u|  at  t=%.0f  (Re=%.0f, %d×%d grid)",
                                p.T, p.Re, N, N),
                 colorbar_title="|u|")
    contour!(p2, xg, yg, ψ', levels=18, color=:white, lw=0.7, alpha=0.85)

    plt = plot(p1, p2, layout=(1,2), size=size, dpi=dpi,
               plot_title="Navier-Stokes:  ∂u/∂t + (u·∇)u = −∇p/ρ + ν∇²u   (incompressible)",
               plot_titlefontsize=13)
    isempty(savepath) || savefig(plt, savepath)
    return plt
end

# ══════════════════════════════════════════════════════════════════════════════
# Generic dispatch on PDESolution
# ══════════════════════════════════════════════════════════════════════════════

"""
    plot_solution(sol::PDESolution; kwargs...) -> Plot

Produce the standard figure for the problem type embedded in `sol.problem`.

| `sol.problem` type       | Figure produced                               |
|--------------------------|-----------------------------------------------|
| `HeatEquation{1}`        | 2-panel: snapshots + space-time heatmap       |
| `WaveEquation{1}`        | 2-panel: snapshots + space-time heatmap       |
| `PoissonEquation{2}`     | 3-panel: f(x,y) / u(x,y) / 3-D surface       |
| `LidCavityFlow`          | 2-panel: schematic + speed heatmap            |
| everything else          | simple heatmap (2-D) or line plot (1-D)       |

All keyword arguments are forwarded to the underlying `fig_*` function
(`savepath`, `size`, `dpi`, and problem-specific options).
"""
function plot_solution(sol::PDESolution{T,1}; kwargs...) where T
    p = sol.problem
    p isa HeatEquation && return fig_heat_equation(p; kwargs...)
    p isa WaveEquation && return fig_wave_equation(p; kwargs...)
    return plot(sol.x, sol.u; xlabel="x", ylabel="u",
                title="$(typeof(p))  t=$(sol.t)")
end

function plot_solution(sol::PDESolution{T,2}; kwargs...) where T
    p = sol.problem
    p isa PoissonEquation && return fig_poisson_equation(p; kwargs...)
    p isa LidCavityFlow   && return fig_navier_stokes(p; kwargs...)
    return heatmap(sol.x, sol.y, sol.u'; xlabel="x", ylabel="y",
                   title="$(typeof(p))  t=$(sol.t)")
end

# ══════════════════════════════════════════════════════════════════════════════
# Extend Plots.plot for Problem structs and PDESolution
# — lets users call plot(prob) or plot(sol) directly
# ══════════════════════════════════════════════════════════════════════════════

Plots.plot(p::HeatEquation{1,F};             kwargs...) where F         = fig_heat_equation(p; kwargs...)
Plots.plot(p::WaveEquation{1,F};             kwargs...) where F         = fig_wave_equation(p; kwargs...)
Plots.plot(p::PoissonEquation{2,F,EF};       kwargs...) where {F,EF}    = fig_poisson_equation(p; kwargs...)
Plots.plot(p::LidCavityFlow;                 kwargs...)                  = fig_navier_stokes(p; kwargs...)
Plots.plot(sol::PDESolution{T,1};            kwargs...) where T          = plot_solution(sol; kwargs...)
Plots.plot(sol::PDESolution{T,2};            kwargs...) where T          = plot_solution(sol; kwargs...)
