using Pkg
Pkg.activate(".")
using juliaPDEs
using CairoMakie
using Printf, Statistics

# ── Problem: 1-D heat equation with a separable IC whose exact solution
#    is just one decaying Fourier mode:
#
#        ∂u/∂t = α ∂²u/∂x²,    x ∈ [0,1],   u(0,t) = u(1,t) = 0
#        u(x, 0) = sin(πx)
#
#        ⇒  u(x, t) = sin(πx) · exp(-α π² t)
#
# Verifying that the numerical solution matches this exact answer is
# *the* sanity check that the discretisation is implemented correctly.
const α = 0.05
exact(x, t) = sin(π * x) * exp(-α * π^2 * t)

# Solve to time T on an N-point endpoint grid. We use Crank–Nicolson
# (θ = 0.5, unconditionally stable, O(Δt²) in time) so the temporal
# error stays well below the spatial error at every N — that lets the
# log-log slope in panel B come out cleanly at the expected 2 instead
# of being polluted by the FE temporal floor when h is small but Δt
# is held constant.
const NT_CN = 400          # Δt = T/NT_CN ≈ 0.0025–0.005 → temporal err ≪ spatial

function solve_to_fe(T_end::Float64, N::Int)
    # Forward Euler — only used for the snapshot panel as a visual sanity check.
    h      = 1.0 / N
    dt_cap = 0.4 * h^2 / α
    Nt     = max(500, ceil(Int, T_end / dt_cap))
    grid   = Grid(a = (0.0,), b = (1.0,), stepsize = (h,))
    tg     = TestGrid(grid, ((t,) -> 0.0,), ((t,) -> 0.0,), (x,) -> sin(π * x))
    return solve(HeatEquation(testgrid = tg, Nt = Nt, T = T_end, α = α))
end

function solve_to_cn(T_end::Float64, N::Int; Nt::Int = NT_CN)
    h    = 1.0 / N
    grid = Grid(a = (0.0,), b = (1.0,), stepsize = (h,))
    tg   = TestGrid(grid, ((t,) -> 0.0,), ((t,) -> 0.0,), (x,) -> sin(π * x))
    return solve_implicit(HeatEquation(testgrid = tg, Nt = Nt, T = T_end, α = α);
                          θ = 0.5)
end

# ── Convergence study at fixed t = 1.0 ──────────────────────────────
Ns       = [16, 32, 64, 128, 256]
errors   = Float64[]
hs       = Float64[]
println("Convergence study at t = 1.0  (Crank–Nicolson, Nt = $NT_CN):")
for N in Ns
    sol = solve_to_cn(1.0, N)
    x   = sol.x
    e   = sol.u .- exact.(x, 1.0)
    err = sqrt(sum(e .^ 2) * (x[2] - x[1]))      # discrete L² with width
    push!(errors, err)
    push!(hs, x[2] - x[1])
    @printf("  N = %4d   h = %.4e   ‖u_num − u_exact‖₂ = %.4e\n", N, x[2]-x[1], err)
end

# Empirical convergence order between successive grids
println("\n  empirical orders:")
for k in 2:length(Ns)
    p = log(errors[k-1] / errors[k]) / log(hs[k-1] / hs[k])
    @printf("    N %4d → %4d   order ≈ %.3f\n", Ns[k-1], Ns[k], p)
end

# ── Snapshots at multiple times for the left panel ──────────────────
# Use forward Euler here so the snapshot panel showcases the actual
# §6.1 algorithm.
ts_plot = [0.0, 0.5, 1.0, 1.5, 2.0]
snaps   = Dict{Float64, Any}()
for t in ts_plot
    snaps[t] = t == 0.0 ? nothing : solve_to_fe(t, 200)
end

# ── Build the figure ────────────────────────────────────────────────
set_theme!(theme_minimal())
update_theme!(
    fontsize             = 14,
    backgroundcolor      = :white,
    Axis = (
        backgroundcolor      = :white,
        xgridcolor           = (:gray, 0.18),
        ygridcolor           = (:gray, 0.18),
        xminorgridvisible    = true,
        yminorgridvisible    = true,
        xminorgridcolor      = (:gray, 0.08),
        yminorgridcolor      = (:gray, 0.08),
        spinewidth           = 1.0,
        xtickwidth           = 1.0,
        ytickwidth           = 1.0,
        titlefont            = :regular,
        titlesize            = 16,
        titlegap             = 12,
        xlabelfont           = :regular,
        ylabelfont           = :regular,
        xlabelsize           = 14,
        ylabelsize           = 14,
    ),
)

fig = Figure(size = (1200, 480))

# Panel A — exact vs numerical snapshots ────────────────────────────
ax1 = Axis(fig[1, 1],
    xlabel = "x",
    ylabel = "u(x, t)",
    title  = "Verification — analytical vs Forward-Euler",
    xticks = 0:0.2:1.0,
)
xlims!(ax1, 0, 1)
ylims!(ax1, -0.05, 1.05)

palette = cgrad(:viridis, length(ts_plot); categorical = true)
xs_fine = range(0, 1, length = 401)

# Plot order: exact line first (dashed), then numerical scatter + thin line.
for (i, t) in enumerate(ts_plot)
    col = palette[i]
    # Exact (dashed)
    lines!(ax1, xs_fine, exact.(xs_fine, t),
           color = col, linewidth = 2.4, linestyle = :dash)
    if t == 0.0
        # IC — solid line stand-in for the "numerical" curve (same as exact)
        lines!(ax1, xs_fine, exact.(xs_fine, t),
               color = col, linewidth = 1.5,
               label = @sprintf("t = %.1f", t))
    else
        sol = snaps[t]
        # Numerical: line + sparse markers so the overlay with the exact dash is visible
        scatter!(ax1, sol.x[1:10:end], sol.u[1:10:end],
                 color = col, markersize = 8, marker = :circle, strokewidth = 0)
        lines!(ax1, sol.x, sol.u,
               color = col, linewidth = 1.5,
               label = @sprintf("t = %.1f", t))
    end
end

axislegend(ax1, position = :rt, framevisible = false,
           labelsize = 12, rowgap = 2, patchsize = (18, 12))
text!(ax1, 0.5, -0.015;
      text = "dashed = exact    ·    solid + dots = numerical",
      align = (:center, :bottom), fontsize = 11, color = :gray30)

# Panel B — convergence ─────────────────────────────────────────────
ax2 = Axis(fig[1, 2],
    xlabel = "grid points  N",
    ylabel = "L² error at t = 1",
    xscale = log10, yscale = log10,
    title  = "Convergence — second-order in Δx",
    xticks = (Ns, string.(Ns)),
)

# Measured error
scatterlines!(ax2, Ns, errors,
              color = :crimson, marker = :diamond, markersize = 12,
              linewidth = 2.2, label = "measured")

# O(h²) reference, normalised at the coarsest grid
ref = errors[1] .* (hs ./ hs[1]) .^ 2
lines!(ax2, Ns, ref,
       color = :gray40, linestyle = :dash, linewidth = 1.6,
       label = "O(Δx²) reference")

axislegend(ax2, position = :lb, framevisible = false, labelsize = 12)

# Print measured slope on the panel
order = log(errors[end] / errors[1]) / log(hs[end] / hs[1])
text!(ax2, 30, errors[1] / 6;
      text = @sprintf("slope ≈ %.2f  ✓", order),
      align = (:left, :center),
      fontsize = 13, color = :crimson, font = :bold)

# Slight padding around the data
ylims!(ax2, minimum(errors) / 4, maximum(errors) * 2)

# ── Save ────────────────────────────────────────────────────────────
outdir  = joinpath(@__DIR__, "..", "site", "output", "figures")
mkpath(outdir)
outpath = joinpath(outdir, "heat_verification.png")
save(outpath, fig; px_per_unit = 2)
println("\nSaved $outpath")
