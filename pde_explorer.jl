### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity
# (so this cell adds the dependency on PlutoUI)
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a1000001-0000-0000-0000-000000000001
md"""
# Julia for HPC — Interactive PDE Explorer

A companion notebook to the **Julia for HPC Student Setup Guide**.
Drag the sliders to change parameters and watch the solutions update in real time.

> **Setup**: This notebook requires Pluto.jl. Install it with:
> ```julia
> using Pkg; Pkg.add("Pluto"); Pkg.add("PlutoUI"); Pkg.add("Plots")
> using Pluto; Pluto.run()
> ```
> Then open this `.jl` file from the Pluto interface.
"""

# ╔═╡ a1000002-0000-0000-0000-000000000002
using PlutoUI, Plots, LinearAlgebra, SparseArrays, Printf

# ╔═╡ a1000003-0000-0000-0000-000000000003
md"""
## 🌡️ 1. Heat Equation (Parabolic)

$$\frac{\partial u}{\partial t} = \alpha \nabla^2 u, \qquad u(0,t) = u(L,t) = 0$$

Use the sliders to change the diffusivity, grid resolution, and snapshot time.
"""

# ╔═╡ a1000004-0000-0000-0000-000000000004
md"""
**Diffusivity** α: $(@bind α_heat Slider(0.001:0.001:0.05, default=0.01, show_value=true))

**Grid points** N: $(@bind N_heat Slider(50:25:400, default=200, show_value=true))

**Snapshot time** t: $(@bind t_heat Slider(0.0:0.05:1.5, default=0.2, show_value=true))

**Initial pulse width**: $(@bind pulse_w Slider(20:20:300, default=100, show_value=true))
"""

# ╔═╡ a1000005-0000-0000-0000-000000000005
function solve_heat_at_time(N, α, t_target; L=1.0, pulse_w=100)
    dx = L / (N+1)
    dt = 0.4 * dx^2 / α
    x  = collect(range(dx, L-dx, length=N))
    u  = exp.(-pulse_w .* (x .- 0.5).^2)
    nsteps = max(1, Int(ceil(t_target / dt)))
    for _ in 1:nsteps
        u[2:end-1] .+= dt * α / dx^2 .*
            (u[3:end] .- 2 .* u[2:end-1] .+ u[1:end-2])
    end
    return x, u
end

# ╔═╡ a1000006-0000-0000-0000-000000000006
let
    x, u = solve_heat_at_time(N_heat, α_heat, t_heat; pulse_w=pulse_w)
    _, u0 = solve_heat_at_time(N_heat, α_heat, 0.0; pulse_w=pulse_w)
    plot(x, u0, lw=2, color=:gray, ls=:dash, label="t = 0",
         xlabel="x", ylabel="u(x, t)",
         title=@sprintf("Heat eq.  α=%.3f, N=%d, t=%.2f", α_heat, N_heat, t_heat),
         ylim=(-0.05, 1.05), framestyle=:box, grid=true, gridalpha=0.3)
    plot!(x, u, lw=2.5, color=:darkorange, label=@sprintf("t = %.2f", t_heat))
end

# ╔═╡ a1000007-0000-0000-0000-000000000007
md"""
### Key observations to discuss with students:
- **Diffusivity α controls speed**: higher α → faster smoothing.
- **CFL stability**: dt = 0.4·dx²/α. If you violated this (try the formula at home), the solution would explode.
- **Grid refinement**: increase N — solution should converge to a smooth limit.
"""

# ╔═╡ a2000001-0000-0000-0000-000000000001
md"""
## 🌊 2. Wave Equation (Hyperbolic)

$$\frac{\partial^2 u}{\partial t^2} = c^2 \nabla^2 u, \qquad u(0,t) = u(L,t) = 0$$

Watch how the wave preserves features (unlike the heat equation, which smooths them).
"""

# ╔═╡ a2000002-0000-0000-0000-000000000002
md"""
**Wave speed** c: $(@bind c_wave Slider(0.5:0.1:3.0, default=1.0, show_value=true))

**Grid points** N: $(@bind N_wave Slider(50:50:500, default=300, show_value=true))

**Snapshot time** t: $(@bind t_wave Slider(0.0:0.05:2.0, default=0.5, show_value=true))

**Initial shape**: $(@bind wave_shape Select(["Triangular pluck", "Sinusoid", "Gaussian"]))
"""

# ╔═╡ a2000003-0000-0000-0000-000000000003
function solve_wave_at_time(N, c, t_target; L=1.0, shape="Triangular pluck")
    dx = L / (N+1)
    dt = 0.4 * dx / c
    x  = collect(range(dx, L-dx, length=N))
    λ  = (c*dt/dx)^2

    u_prev = if shape == "Triangular pluck"
        [xi < 0.5 ? 2xi : 2(1-xi) for xi in x]
    elseif shape == "Sinusoid"
        sin.(π .* x)
    else  # Gaussian
        exp.(-100 .* (x .- 0.5).^2)
    end
    u_curr = copy(u_prev)
    u_next = similar(u_curr)

    nsteps = max(0, Int(ceil(t_target / dt)))
    for _ in 1:nsteps
        u_next[2:end-1] .= 2 .* u_curr[2:end-1] .- u_prev[2:end-1] .+
            λ .* (u_curr[3:end] .- 2 .* u_curr[2:end-1] .+ u_curr[1:end-2])
        u_prev, u_curr = u_curr, copy(u_next)
    end
    return x, u_curr, u_prev
end

# ╔═╡ a2000004-0000-0000-0000-000000000004
let
    x, u, _ = solve_wave_at_time(N_wave, c_wave, t_wave; shape=wave_shape)
    _, u0, _ = solve_wave_at_time(N_wave, c_wave, 0.0; shape=wave_shape)
    plot(x, u0, lw=2, color=:gray, ls=:dash, label="t = 0",
         xlabel="x", ylabel="u(x, t)",
         title=@sprintf("Wave eq.  c=%.1f, N=%d, t=%.2f  (%s)",
                        c_wave, N_wave, t_wave, wave_shape),
         ylim=(-1.2, 1.2), framestyle=:box, grid=true, gridalpha=0.3)
    plot!(x, u, lw=2.5, color=:steelblue, label=@sprintf("t = %.2f", t_wave))
    hline!([0], color=:gray, lw=0.5, label="")
end

# ╔═╡ a2000005-0000-0000-0000-000000000005
md"""
### Key observations:
- **Wave splits**: a triangular pluck splits into two waves moving in opposite directions.
- **Reflections**: each wave bounces off the boundary at x=0 and x=L.
- **No dissipation**: the wave equation conserves energy — try long times!
- **CFL**: dt < dx/c. Less restrictive than the heat equation.
"""

# ╔═╡ a3000001-0000-0000-0000-000000000001
md"""
## ⚡ 3. Poisson Equation (Elliptic)

$$-\nabla^2 u = f, \qquad u = 0 \text{ on the boundary}$$

No time evolution — pick a source f, get a steady-state u in one solve.
"""

# ╔═╡ a3000002-0000-0000-0000-000000000002
md"""
**Grid points per side** N: $(@bind N_poi Slider(20:10:80, default=50, show_value=true))

**Source mode** (m,n): $(@bind mode_m Slider(1:5, default=1, show_value=true)), $(@bind mode_n Slider(1:5, default=1, show_value=true))
"""

# ╔═╡ a3000003-0000-0000-0000-000000000003
function solve_poisson_modes(N, m, n; L=1.0)
    h  = L / (N+1)
    xs = collect(range(h, L-h, length=N))
    ys = collect(range(h, L-h, length=N))
    e   = ones(N)
    T1D = spdiagm(0 => -2 .* e, 1 => e[1:end-1], -1 => e[1:end-1])
    I_N = sparse(I, N, N)
    A   = -(kron(I_N, T1D) + kron(T1D, I_N)) ./ h^2
    F = [(m^2 + n^2)*π^2 * sin(m*π*x) * sin(n*π*y) for y in ys, x in xs]
    u = A \ vec(F)
    return xs, ys, reshape(u, N, N), F
end

# ╔═╡ a3000004-0000-0000-0000-000000000004
let
    xs, ys, U, F = solve_poisson_modes(N_poi, mode_m, mode_n)
    p1 = heatmap(xs, ys, F, c=:coolwarm, aspect_ratio=:equal,
                 title="Source f(x,y)", xlabel="x", ylabel="y")
    p2 = heatmap(xs, ys, U, c=:viridis, aspect_ratio=:equal,
                 title="Solution u(x,y)", xlabel="x", ylabel="y")
    plot(p1, p2, layout=(1,2), size=(800, 350),
         plot_title=@sprintf("Poisson modes (m,n) = (%d,%d), N=%d",
                             mode_m, mode_n, N_poi),
         plot_titlefontsize=12)
end

# ╔═╡ a3000005-0000-0000-0000-000000000005
md"""
### Key observations:
- **Mode (1,1)** gives a single bump — the lowest eigenmode.
- **Higher modes** give grid patterns of (m × n) cells.
- **No time** — this is solved in a single sparse linear solve `u = A \\ b`.
- The condition number of A grows like 1/h², so very fine grids may need iterative solvers + preconditioning.
"""

# ╔═╡ a4000001-0000-0000-0000-000000000001
md"""
## 📈 4. Convergence Study (Live)

Watch the L∞ error decrease as we refine the grid. The **slope on a log-log plot reveals the order of accuracy**.
"""

# ╔═╡ a4000002-0000-0000-0000-000000000002
md"""
**Maximum N**: $(@bind N_max_conv Slider([20, 40, 80, 160, 320], default=160, show_value=true))
"""

# ╔═╡ a4000003-0000-0000-0000-000000000003
let
    Ns_all = [10, 20, 40, 80, 160, 320]
    Ns = filter(n -> n <= N_max_conv, Ns_all)
    hs, errs, ts = Float64[], Float64[], Float64[]
    for N in Ns
        h = 1.0/(N+1)
        t0 = time()
        xs, ys, U, _ = solve_poisson_modes(N, 1, 1)
        t1 = time()
        Ue = [sin(π*x)*sin(π*y) for y in ys, x in xs]
        push!(hs, h); push!(errs, maximum(abs.(U .- Ue)))
        push!(ts, t1-t0)
    end
    ref = errs[1] .* (hs ./ hs[1]).^2

    p1 = plot(hs, errs, marker=:circle, ms=7, lw=2.2,
              xscale=:log10, yscale=:log10,
              label="L∞ error",
              xlabel="grid spacing h", ylabel="error",
              title="Convergence", framestyle=:box, grid=true, gridalpha=0.3)
    plot!(p1, hs, ref, ls=:dash, lw=1.5, label="O(h²) reference", color=:firebrick)

    p2 = plot(Ns, ts, marker=:diamond, ms=7, lw=2.2, color=:seagreen,
              xscale=:log10, yscale=:log10,
              xlabel="N", ylabel="solve time (s)",
              title="Runtime", framestyle=:box, grid=true, gridalpha=0.3,
              legend=false)

    plot(p1, p2, layout=(1,2), size=(800, 350))
end

# ╔═╡ a4000004-0000-0000-0000-000000000004
md"""
### Reading the convergence plot:
- The numerical curve should **lie parallel to the O(h²) reference** for a 2nd-order scheme.
- If your slope is **less steep**, your code has a bug or is reverting to 1st-order at the boundary.
- If your slope is **steeper**, you may have stumbled into superconvergence (rare, usually requires special grids).

This is the most important diagnostic in numerical PDE work — always run a convergence study!
"""

# ╔═╡ a5000001-0000-0000-0000-000000000001
md"""
## 🎓 Discussion Prompts for Class

1. **Heat eq.** — what happens if you set α much larger? Why does the simulation become unstable? (Answer: violates CFL, dt becomes too large.)
2. **Wave eq.** — increase wave speed c and watch the CFL margin shrink. Try c above the safe limit.
3. **Poisson** — pick mode (3, 2) and discuss why the solution looks like the source but smoother.
4. **Convergence** — what happens if you accidentally use a 1st-order finite difference at one boundary? Predict the slope.

---

> 📚 For the full setup guide and HPC/Slurm instructions, see `julia_hpc_setup_guide.docx`.
> 🔧 To regenerate all figures and the docx, run: `julia --project=. generate_lecture_doc.jl`
"""

# ╔═╡ Cell order:
# ╟─a1000001-0000-0000-0000-000000000001
# ╠═a1000002-0000-0000-0000-000000000002
# ╟─a1000003-0000-0000-0000-000000000003
# ╟─a1000004-0000-0000-0000-000000000004
# ╠═a1000005-0000-0000-0000-000000000005
# ╠═a1000006-0000-0000-0000-000000000006
# ╟─a1000007-0000-0000-0000-000000000007
# ╟─a2000001-0000-0000-0000-000000000001
# ╟─a2000002-0000-0000-0000-000000000002
# ╠═a2000003-0000-0000-0000-000000000003
# ╠═a2000004-0000-0000-0000-000000000004
# ╟─a2000005-0000-0000-0000-000000000005
# ╟─a3000001-0000-0000-0000-000000000001
# ╟─a3000002-0000-0000-0000-000000000002
# ╠═a3000003-0000-0000-0000-000000000003
# ╠═a3000004-0000-0000-0000-000000000004
# ╟─a3000005-0000-0000-0000-000000000005
# ╟─a4000001-0000-0000-0000-000000000001
# ╟─a4000002-0000-0000-0000-000000000002
# ╠═a4000003-0000-0000-0000-000000000003
# ╟─a4000004-0000-0000-0000-000000000004
# ╟─a5000001-0000-0000-0000-000000000001
