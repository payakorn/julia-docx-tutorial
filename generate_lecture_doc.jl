#!/usr/bin/env julia
# ==============================================================================
#  generate_lecture_doc.jl
#
#  Generates the "Julia for HPC — Student Setup Guide" document, including
#  PDE figures (Heat, Wave, Poisson, Navier-Stokes) and a full .docx file.
#
#  Usage:
#     julia --project=. generate_lecture_doc.jl
#
#  All knobs are in the CONFIG block below — edit those, re-run, done.
#
#  Required packages:
#     Plots, LinearAlgebra, SparseArrays, Printf
#  Optional (auto-detected): PythonCall (preferred) or PyCall, for python-docx
# ==============================================================================

using Pkg
# ---------------- Auto-install missing packages on first run ------------------
const REQUIRED_PKGS = ["Plots", "PythonCall", "Gmsh", "CairoMakie"]

function ensure_packages()
    installed = Set(keys(Pkg.project().dependencies))
    for p in REQUIRED_PKGS
        if !(p in installed)
            @info "Installing missing package: $p"
            Pkg.add(p)
        end
    end
end
ensure_packages()

using Plots, LinearAlgebra, SparseArrays, Printf, Dates
using PythonCall
using Gmsh: gmsh
# Makie is imported under an alias so its exports (heatmap, surface, contour,
# contourf, plot, …) never collide with the identically-named Plots functions
# used elsewhere in this file. Every Makie call below is qualified with `CM.`.
import CairoMakie as CM

# Set headless backend for HPC / CI
ENV["GKSwstype"] = "100"
gr()

# =============================================================================
# CONFIG  — Edit this block to change any parameter
# =============================================================================
Base.@kwdef mutable struct Config
    # ---- output paths ----
    output_dir   :: String = "output"            # where the .docx lands
    figures_dir  :: String = "site/output/figures"  # served at /output/figures/ on the site
    docx_name    :: String = "julia_hpc_setup_guide.docx"

    # ---- figure styling ----
    fig_dpi      :: Int    = 150
    fig_width    :: Int    = 1100   # px
    fig_height   :: Int    = 450    # px

    # ---- HEAT EQUATION ----
    heat_N         :: Int     = 200
    heat_alpha     :: Float64 = 0.01
    heat_L         :: Float64 = 1.0
    heat_T_final   :: Float64 = 1.0
    heat_cfl_safety:: Float64 = 0.4         # dt = safety * dx² / α
    heat_snapshots :: Vector{Float64} = [0.0, 0.05, 0.2, 0.5, 1.0]
    heat_init      :: Function = x -> exp(-100*(x-0.5)^2)   # Gaussian pulse

    # ---- WAVE EQUATION ----
    wave_N         :: Int     = 300
    wave_c         :: Float64 = 1.0
    wave_L         :: Float64 = 1.0
    wave_T_final   :: Float64 = 1.5
    wave_T_show    :: Float64 = 2.0          # for space-time plot
    wave_cfl_safety:: Float64 = 0.4
    wave_snapshots :: Vector{Float64} = [0.1, 0.3, 0.5, 0.8, 1.2]
    wave_init      :: Function = x -> x < 0.5 ? 2x : 2*(1-x)  # triangular pluck

    # ---- POISSON EQUATION ----
    poisson_N      :: Int     = 50
    poisson_L      :: Float64 = 1.0
    poisson_f      :: Function = (x,y) -> 2*π^2 * sin(π*x) * sin(π*y)
    poisson_exact  :: Function = (x,y) -> sin(π*x) * sin(π*y)

    # ---- NAVIER-STOKES (lid-driven cavity sketch) ----
    ns_Re          :: Float64 = 100.0
    ns_U_lid       :: Float64 = 1.0
    ns_grid_n      :: Int     = 40    # for streamline visualization

    # ---- doc title-page text ----
    doc_title      :: String = "Julia for HPC"
    doc_subtitle   :: String = "Student Setup Guide"
    doc_tagline    :: String = "WSL • juliaup • Julia • Packages • PDE Solvers • Slurm"
    doc_author     :: String = "Course Instructor"
end

const CFG = Config()

# Make output dirs
mkpath(CFG.output_dir)
mkpath(CFG.figures_dir)

# Color palette
const BLUE    = "#1F4E79"
const ACCENT  = "#2E75B6"
const ORANGE  = "#E07B39"
const GREEN   = "#588157"
const RED     = "#C84B31"

# =============================================================================
# PDE SOLVERS  — pure Julia, used both in the lecture and for figure generation
# =============================================================================

"""
    solve_heat_1d(; N, α, L, T_final, init_func, cfl_safety, snapshot_times)

Solves ∂u/∂t = α ∇²u on [0, L] with u(0)=u(L)=0 (Dirichlet).
Returns (x, snapshots_dict, full_history, t_grid).
"""
function solve_heat_1d(; N=CFG.heat_N, α=CFG.heat_alpha, L=CFG.heat_L,
                       T_final=CFG.heat_T_final, init_func=CFG.heat_init,
                       cfl_safety=CFG.heat_cfl_safety,
                       snapshot_times=CFG.heat_snapshots)
    dx = L / (N+1)
    dt = cfl_safety * dx^2 / α
    x  = collect(range(dx, L-dx, length=N))
    nsteps = Int(ceil(T_final / dt))

    u = init_func.(x)
    snapshots = Dict{Float64, Vector{Float64}}()
    snapshots[0.0] = copy(u)

    # full history for space-time plot
    full_hist = [copy(u)]
    t_grid    = [0.0]
    save_every = max(1, nsteps ÷ 200)

    target_times = sort(filter(t -> t > 0, snapshot_times))
    t = 0.0
    for step in 1:nsteps
        u[2:end-1] .+= dt * α / dx^2 .*
            (u[3:end] .- 2 .* u[2:end-1] .+ u[1:end-2])
        t += dt
        for tt in target_times
            if !haskey(snapshots, tt) && t >= tt
                snapshots[tt] = copy(u)
            end
        end
        if step % save_every == 0
            push!(full_hist, copy(u))
            push!(t_grid, t)
        end
    end
    return x, snapshots, hcat(full_hist...), t_grid
end

"""
    solve_wave_1d(; N, c, L, T_final, init_func, cfl_safety, snapshot_times)

Solves ∂²u/∂t² = c² ∇²u on [0, L] with u(0)=u(L)=0 and zero initial velocity.
Uses leap-frog scheme.
"""
function solve_wave_1d(; N=CFG.wave_N, c=CFG.wave_c, L=CFG.wave_L,
                       T_final=CFG.wave_T_final, init_func=CFG.wave_init,
                       cfl_safety=CFG.wave_cfl_safety,
                       snapshot_times=CFG.wave_snapshots)
    dx = L / (N+1)
    dt = cfl_safety * dx / c
    x  = collect(range(dx, L-dx, length=N))
    λ  = (c*dt/dx)^2
    nsteps = Int(ceil(T_final / dt))

    u_prev = init_func.(x)
    u_curr = copy(u_prev)        # zero initial velocity
    u_next = similar(u_curr)

    snapshots = Dict{Float64, Vector{Float64}}()
    snapshots[0.0] = copy(u_curr)

    full_hist = [copy(u_curr)]
    t_grid    = [0.0]
    save_every = max(1, nsteps ÷ 250)

    target_times = sort(filter(t -> t > 0, snapshot_times))
    t = 0.0
    for step in 1:nsteps
        u_next[2:end-1] .= 2 .* u_curr[2:end-1] .- u_prev[2:end-1] .+
            λ .* (u_curr[3:end] .- 2 .* u_curr[2:end-1] .+ u_curr[1:end-2])
        u_prev, u_curr = u_curr, copy(u_next)
        t += dt
        for tt in target_times
            if !haskey(snapshots, tt) && t >= tt
                snapshots[tt] = copy(u_curr)
            end
        end
        if step % save_every == 0
            push!(full_hist, copy(u_curr))
            push!(t_grid, t)
        end
    end
    return x, snapshots, hcat(full_hist...), t_grid
end

"""
    solve_poisson_2d(; N, L, f_func)

Solves -∇²u = f on [0,L]² with u=0 on the boundary, using a 5-point
finite-difference stencil and sparse direct solve.
"""
function solve_poisson_2d(; N=CFG.poisson_N, L=CFG.poisson_L,
                          f_func=CFG.poisson_f)
    h  = L / (N+1)
    xs = collect(range(h, L-h, length=N))
    ys = collect(range(h, L-h, length=N))

    # Build 2D Laplacian via Kronecker products
    e   = ones(N)
    T1D = spdiagm(0 => -2 .* e, 1 => e[1:end-1], -1 => e[1:end-1])
    I_N = sparse(I, N, N)
    A   = -(kron(I_N, T1D) + kron(T1D, I_N)) ./ h^2

    F = [f_func(x, y) for y in ys, x in xs]
    b = vec(F)

    u = A \ b
    U = reshape(u, N, N)

    return xs, ys, U, F
end

# =============================================================================
# FIGURE GENERATION  — produces all 6 figures used in the document
# =============================================================================

function fig_pde_pipeline()
    @info "Generating Figure 1: PDE pipeline diagram"
    stages = [("Model",      "PDE + BCs + ICs\nClassify type"),
              ("Discretize", "Grid, FDM/FEM\nTime scheme"),
              ("Build System","Assemble A, b\nSparse matrices"),
              ("Solve",      "Direct / iterative\nTime-stepping"),
              ("Visualize",  "Plot, verify\nConvergence")]
    colors = [:steelblue, :dodgerblue, :seagreen, :darkorange, :gray]

    plt = plot(size=(1200, 380), legend=false, framestyle=:none,
               xlim=(-0.3, 12.5), ylim=(0, 2.2), aspect_ratio=:equal,
               background_color=:white, dpi=CFG.fig_dpi)

    for (i, (title, desc)) in enumerate(stages)
        x = (i-1) * 2.4
        # box
        plot!(plt, [x, x+2, x+2, x, x], [0.3, 0.3, 1.7, 1.7, 0.3],
              seriestype=:shape, fillcolor=colors[i], linecolor=:black,
              linewidth=1.5, fillalpha=0.85)
        annotate!(plt, x+1.0, 1.35,
                  text("$i. $title", :white, :center, 13, "Computer Modern", :bold))
        annotate!(plt, x+1.0, 0.78,
                  text(desc, :white, :center, 10, "Computer Modern"))
        if i < length(stages)
            annotate!(plt, x+2.2, 1.0, text("→", :black, :center, 22, :bold))
        end
    end
    title!(plt, "The PDE-Solving Pipeline")
    savefig(plt, joinpath(CFG.figures_dir, "pde_pipeline.png"))
end

function fig_pde_types()
    @info "Generating Figure 2: PDE type classification"
    types = [("Parabolic",   "Heat / Diffusion\n∂u/∂t = α∇²u",
              "Smooths features\nover time", :darkorange),
             ("Hyperbolic",  "Wave / Advection\n∂²u/∂t² = c²∇²u",
              "Preserves features\nfinite speed c", :dodgerblue),
             ("Elliptic",    "Poisson / Laplace\n−∇²u = f",
              "Steady state\nno time evolution", :seagreen),
             ("Mixed/NL",    "Navier-Stokes\n∂u/∂t + (u·∇)u\n= −∇p + ν∇²u",
              "Hardest!\nFluid dynamics", :firebrick)]

    plt = plot(size=(1200, 500), legend=false, framestyle=:none,
               xlim=(-0.3, 12.0), ylim=(0.8, 4.4),
               background_color=:white, dpi=CFG.fig_dpi)

    for (i, (cat, eq, desc, color)) in enumerate(types)
        x = (i-1) * 3
        plot!(plt, [x, x+2.6, x+2.6, x, x], [1.5, 1.5, 4.0, 4.0, 1.5],
              seriestype=:shape, fillcolor=color, fillalpha=0.18,
              linecolor=color, linewidth=2)
        annotate!(plt, x+1.3, 3.7,
                  text(cat, color, :center, 14, :bold))
        annotate!(plt, x+1.3, 3.0, text(eq, :black, :center, 10))
        annotate!(plt, x+1.3, 2.0,
                  text(desc, RGB(0.35,0.35,0.35), :center, 10, :italic))
    end
    title!(plt, "PDE Classification — Four Major Types")
    savefig(plt, joinpath(CFG.figures_dir, "pde_types.png"))
end

function fig_heat_equation()
    @info "Generating Figure 3: Heat equation"
    x, snapshots, U_hist, t_grid =
        solve_heat_1d()

    times_sorted = sort(collect(keys(snapshots)))
    n = length(times_sorted)
    cmap = cgrad(:YlOrRd, range(0.3, 0.95, length=n))

    p1 = plot(xlabel="x", ylabel="u(x, t)",
              title="Heat Equation:  ∂u/∂t = α ∇²u",
              titlefontcolor=parse(Colorant, BLUE), legend=:topright,
              grid=true, gridalpha=0.3, ylim=(-0.05, 1.05), framestyle=:box)
    for (i, tt) in enumerate(times_sorted)
        plot!(p1, x, snapshots[tt], lw=2.2, color=cmap[i],
              label=@sprintf("t = %.2f", tt))
    end
    # diffusion arrow
    plot!(p1, [0.5, 0.5], [0.95, 0.30], arrow=true, color=:firebrick,
          lw=1.5, label="")
    annotate!(p1, 0.62, 0.65,
              text("diffuses\n& smooths", :firebrick, :left, 9))

    p2 = heatmap(t_grid, x, U_hist, c=:YlOrRd,
                 xlabel="time t", ylabel="x",
                 title="Space-time evolution",
                 titlefontcolor=parse(Colorant, BLUE), colorbar_title="u")

    plt = plot(p1, p2, layout=(1,2), size=(CFG.fig_width, CFG.fig_height),
               dpi=CFG.fig_dpi, plot_titlefontsize=13)
    savefig(plt, joinpath(CFG.figures_dir, "heat_equation.png"))
end

function fig_wave_equation()
    @info "Generating Figure 4: Wave equation"
    x, snapshots, U_hist, t_grid = solve_wave_1d()

    times_sorted = sort(collect(keys(snapshots)))
    n = length(times_sorted)
    cmap = cgrad(:viridis, range(0.1, 0.9, length=n))

    p1 = plot(xlabel="x", ylabel="u(x, t)",
              title="Wave Equation:  ∂²u/∂t² = c² ∇²u",
              titlefontcolor=parse(Colorant, BLUE), legend=:topright,
              grid=true, gridalpha=0.3, ylim=(-1.1, 1.1), framestyle=:box)
    hline!(p1, [0], color=:gray, lw=0.5, label="")
    for (i, tt) in enumerate(times_sorted)
        plot!(p1, x, snapshots[tt], lw=2.0, color=cmap[i],
              label=@sprintf("t = %.2f", tt))
    end

    # Run a longer simulation for the space-time plot
    _, _, U_long, t_long = solve_wave_1d(T_final=CFG.wave_T_show)

    p2 = heatmap(t_long, x, U_long, c=:RdBu,
                 clims=(-1, 1),
                 xlabel="time t", ylabel="x",
                 title="Space-time: characteristic lines",
                 titlefontcolor=parse(Colorant, BLUE), colorbar_title="u")

    plt = plot(p1, p2, layout=(1,2), size=(CFG.fig_width, CFG.fig_height),
               dpi=CFG.fig_dpi)
    savefig(plt, joinpath(CFG.figures_dir, "wave_equation.png"))
end

function fig_poisson_equation()
    @info "Generating Figure 5: Poisson equation"
    xs, ys, U, F = solve_poisson_2d()

    p1 = heatmap(xs, ys, F, c=:coolwarm,
                 xlabel="x", ylabel="y",
                 title="Source term  f(x,y)",
                 titlefontcolor=parse(Colorant, BLUE),
                 aspect_ratio=:equal)
    p2 = heatmap(xs, ys, U, c=:viridis,
                 xlabel="x", ylabel="y",
                 title="Solution  u(x,y)",
                 titlefontcolor=parse(Colorant, BLUE),
                 aspect_ratio=:equal)
    p3 = surface(xs, ys, U, c=:viridis, camera=(40, 30),
                 xlabel="x", ylabel="y", zlabel="u",
                 title="3D view of u(x,y)",
                 titlefontcolor=parse(Colorant, BLUE))

    plt = plot(p1, p2, p3, layout=(1,3),
               size=(1300, 480), dpi=CFG.fig_dpi,
               plot_title="Poisson Equation:  −∇²u = f   (steady-state, no time)",
               plot_titlefontsize=14,
               plot_titlefontcolor=parse(Colorant, BLUE))
    savefig(plt, joinpath(CFG.figures_dir, "poisson_equation.png"))
end

"""
    solve_lid_cavity_psi_omega(; N, Re, T, U_lid)

Lid-driven cavity at moderate Reynolds number using the streamfunction-
vorticity formulation. Semi-implicit time stepping (implicit diffusion,
explicit convection) so the only stability constraint is the
convective CFL `U·dt/h ≤ 1`. Returns N×N arrays where `arr[i, j]` is
the value at grid point `(x_i, y_j)`, plus the grid spacing `h`.

Equations:
    ω = -∇²ψ                                 (vorticity from streamfunction)
    ∂ω/∂t + (∂ψ/∂y)(∂ω/∂x) - (∂ψ/∂x)(∂ω/∂y) = ν∇²ω    (vorticity transport)
    u =  ∂ψ/∂y,   v = -∂ψ/∂x                (velocity from ψ)

Boundary conditions:
    ψ = 0 on every wall  (no normal flow)
    Thom's formula gives the wall vorticity from ψ at the first interior layer
    Top wall additionally has  u = U_lid  (moving lid)
"""
function solve_lid_cavity_psi_omega(; N::Int=48, Re::Float64=100.0,
                                    T::Float64=10.0, U_lid::Float64=1.0)
    L = 1.0
    h = L / (N - 1)
    ν = U_lid * L / Re

    # Convection-limited time step (diffusion is implicit).
    dt     = 0.30 * h / U_lid
    nsteps = ceil(Int, T / dt)

    ψ = zeros(N, N)
    ω = zeros(N, N)

    # Build the discrete 2-D Laplacian on the (N-2)×(N-2) interior, with
    # Dirichlet zero implicitly assumed at the boundaries.
    n_int = N - 2
    nn    = n_int^2
    II = Int[];  JJ = Int[];  VV = Float64[]
    for jj in 1:n_int, ii in 1:n_int
        k = (jj - 1) * n_int + ii
        push!(II, k); push!(JJ, k); push!(VV, -4.0)
        if ii > 1     ; push!(II, k); push!(JJ, k - 1);     push!(VV, 1.0); end
        if ii < n_int ; push!(II, k); push!(JJ, k + 1);     push!(VV, 1.0); end
        if jj > 1     ; push!(II, k); push!(JJ, k - n_int); push!(VV, 1.0); end
        if jj < n_int ; push!(II, k); push!(JJ, k + n_int); push!(VV, 1.0); end
    end
    Lap = sparse(II, JJ, VV, nn, nn) ./ h^2

    Inn    = sparse(I, nn, nn)
    A_diff = lu(Inn - dt * ν * Lap)        # for vorticity diffusion
    Lap_lu = lu(Lap)                       # for streamfunction Poisson
    coeff  = ν * dt / h^2                  # near-wall stencil weight

    for _ in 1:nsteps
        # 1. Wall vorticity from ψ via Thom's formula
        @views ω[2:N-1, N] .= -2 .* ψ[2:N-1, N-1] ./ h^2 .- 2 * U_lid / h  # top lid
        @views ω[2:N-1, 1] .= -2 .* ψ[2:N-1, 2]   ./ h^2                    # bottom
        @views ω[1, 2:N-1] .= -2 .* ψ[2,   2:N-1] ./ h^2                    # left
        @views ω[N, 2:N-1] .= -2 .* ψ[N-1, 2:N-1] ./ h^2                    # right

        # 2. Explicit advection at every interior node
        ψy = (ψ[2:N-1, 3:N]   .- ψ[2:N-1, 1:N-2]) ./ (2h)
        ψx = (ψ[3:N,   2:N-1] .- ψ[1:N-2, 2:N-1]) ./ (2h)
        ωx = (ω[3:N,   2:N-1] .- ω[1:N-2, 2:N-1]) ./ (2h)
        ωy = (ω[2:N-1, 3:N]   .- ω[2:N-1, 1:N-2]) ./ (2h)
        adv = ψy .* ωx .- ψx .* ωy

        # 3. Build the right-hand side for the implicit diffusion solve
        rhs_int = ω[2:N-1, 2:N-1] .- dt .* adv
        rhs_int[1,   :] .+= coeff .* ω[1,   2:N-1]
        rhs_int[end, :] .+= coeff .* ω[N,   2:N-1]
        rhs_int[:,   1] .+= coeff .* ω[2:N-1, 1]
        rhs_int[:, end] .+= coeff .* ω[2:N-1, N]

        # 4. (I - νΔt L) ω_int^{n+1} = rhs
        ω_int_new = A_diff \ vec(rhs_int)
        ω[2:N-1, 2:N-1] .= reshape(ω_int_new, n_int, n_int)

        # 5. Solve  L ψ = -ω  on the interior  (ψ = 0 on walls)
        ψ_int_new = Lap_lu \ -vec(ω[2:N-1, 2:N-1])
        ψ[2:N-1, 2:N-1] .= reshape(ψ_int_new, n_int, n_int)
    end

    # Recover velocity field by central differences of ψ
    u = zeros(N, N)
    v = zeros(N, N)
    u[:, 2:N-1] .=  (ψ[:, 3:N]   .- ψ[:, 1:N-2])   ./ (2h)
    v[2:N-1, :] .= -(ψ[3:N, :]   .- ψ[1:N-2, :])   ./ (2h)
    u[:, N] .= U_lid                                 # lid velocity

    speed = sqrt.(u.^2 .+ v.^2)
    return ψ, ω, u, v, speed, h
end

function fig_navier_stokes()
    @info "Generating Figure 6: Navier-Stokes (lid-driven cavity, ψ-ω solver)"

    # ---- Left panel: cavity schematic (unchanged) ----
    p1 = plot(legend=false, framestyle=:none, aspect_ratio=:equal,
              xlim=(-0.2, 1.2), ylim=(-0.2, 1.4),
              title="Lid-Driven Cavity (NS benchmark)",
              titlefontcolor=parse(Colorant, BLUE))
    plot!(p1, [0,1,1,0,0], [0,0,1,1,0], color=:black, lw=2.5)
    plot!(p1, [0.15, 0.85], [1.05, 1.05], arrow=true, color=:firebrick, lw=2.5)
    annotate!(p1, 0.5, 1.18,
              text("moving lid  U = $(CFG.ns_U_lid)", :firebrick, :center, 11, :bold))
    annotate!(p1, -0.07, 0.5, text("u=v=0", parse(Colorant, BLUE), :right,  10))
    annotate!(p1,  1.07, 0.5, text("u=v=0", parse(Colorant, BLUE), :left,   10))
    annotate!(p1,  0.5, -0.10, text("u=v=0", parse(Colorant, BLUE), :center, 10))
    θ = range(0, 2π, length=100)
    for r in (0.15, 0.25, 0.35)
        plot!(p1, 0.5 .+ r .* cos.(θ), 0.55 .+ r .* sin.(θ),
              color=:seagreen, lw=1, alpha=0.6)
    end
    annotate!(p1, 0.5, 0.55, text("vortex", :seagreen, :center, 11, :bold))

    # ---- Right panel: REAL numerical solution (ψ-ω solver) ----
    N_grid = 48
    ψ, ω, u, v, speed, h = solve_lid_cavity_psi_omega(
        N=N_grid, Re=CFG.ns_Re, T=10.0, U_lid=CFG.ns_U_lid)
    xg = range(0, 1, length=N_grid)
    yg = range(0, 1, length=N_grid)

    # Heatmap of speed (transposed because Plots' heatmap expects M[y, x]).
    p2 = heatmap(xg, yg, speed', c=:plasma, aspect_ratio=:equal,
                 xlims=(0, 1), ylims=(0, 1),
                 xlabel="x", ylabel="y",
                 title=@sprintf("Numerical |u|  at  t=10  (Re=%.0f, %d×%d grid)",
                                CFG.ns_Re, N_grid, N_grid),
                 titlefontcolor=parse(Colorant, BLUE), colorbar_title="|u|")
    # Overlay streamlines as contours of ψ
    contour!(p2, xg, yg, ψ', levels=18, color=:white, lw=0.7, alpha=0.85)

    plt = plot(p1, p2, layout=(1, 2),
               size=(CFG.fig_width, CFG.fig_height + 50), dpi=CFG.fig_dpi,
               plot_title="Navier-Stokes:  ∂u/∂t + (u·∇)u = −∇p/ρ + ν∇²u   (incompressible)",
               plot_titlefontsize=13,
               plot_titlefontcolor=parse(Colorant, BLUE))
    savefig(plt, joinpath(CFG.figures_dir, "navier_stokes.png"))
end

function fig_convergence_study()
    @info "Generating Figure 7: Convergence study (Poisson)"

    # Run Poisson at multiple resolutions, measure max error vs exact
    Ns = [10, 20, 40, 80, 160]
    hs       = Float64[]
    errors   = Float64[]
    runtimes = Float64[]

    for N in Ns
        h = CFG.poisson_L / (N+1)
        t0 = time()
        xs, ys, U, _ = solve_poisson_2d(N=N)
        t1 = time()
        U_exact = [CFG.poisson_exact(x, y) for y in ys, x in xs]
        push!(hs, h)
        push!(errors, maximum(abs.(U .- U_exact)))
        push!(runtimes, t1 - t0)
    end

    # Reference O(h²) line through the first data point
    ref_h2 = errors[1] .* (hs ./ hs[1]).^2

    p1 = plot(hs, errors, xscale=:log10, yscale=:log10,
              marker=:circle, ms=7, lw=2.2, color=parse(Colorant, ACCENT),
              label="L∞ error (numerical)",
              xlabel="grid spacing h", ylabel="max |u_num − u_exact|",
              title="Convergence — Poisson 2D",
              titlefontcolor=parse(Colorant, BLUE),
              grid=true, gridalpha=0.3, framestyle=:box, legend=:bottomright)
    plot!(p1, hs, ref_h2, ls=:dash, lw=1.8, color=:firebrick,
          label="O(h²) reference")

    # Also plot wall-clock runtime vs N
    p2 = plot(Ns, runtimes, marker=:diamond, ms=7, lw=2.2,
              color=parse(Colorant, GREEN),
              xscale=:log10, yscale=:log10,
              xlabel="grid points per side N", ylabel="solve time (s)",
              title="Sparse direct solver scaling",
              titlefontcolor=parse(Colorant, BLUE),
              grid=true, gridalpha=0.3, framestyle=:box, legend=false)

    plt = plot(p1, p2, layout=(1,2), size=(CFG.fig_width, CFG.fig_height),
               dpi=CFG.fig_dpi)
    savefig(plt, joinpath(CFG.figures_dir, "convergence_study.png"))
end

"""
    gmsh_unit_square(; lc=0.5)

Mesh the unit square [0,1]² with Gmsh and return `(xy, tri_nodes)` where
`xy` is `2 × n_nodes` of node coordinates and `tri_nodes` is `3 × n_tri`
of 1-based global node indices per triangle.
"""
function gmsh_unit_square(; lc::Float64=0.5)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)   # quiet output
    gmsh.model.add("unit_square")

    # 1. Geometry — corner points of [0,1] × [0,1] with target element size lc
    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)
    p2 = gmsh.model.geo.addPoint(1.0, 0.0, 0.0, lc)
    p3 = gmsh.model.geo.addPoint(1.0, 1.0, 0.0, lc)
    p4 = gmsh.model.geo.addPoint(0.0, 1.0, 0.0, lc)

    # 2. Edges and the surface they bound
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    gmsh.model.geo.addPlaneSurface([cl])

    # 3. Synchronise the CAD kernel and triangulate the surface
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)        # 2D ⇒ triangles

    # 4. Pull the mesh back into Julia
    node_tags, coord, _   = gmsh.model.mesh.getNodes()
    _, _, elem_node_tags  = gmsh.model.mesh.getElements(2)

    n_nodes = length(node_tags)
    xy      = reshape(coord, 3, n_nodes)[1:2, :]      # drop z

    # Gmsh tags can be sparse — remap them to dense 1..n_nodes indices.
    tag2idx   = Dict(Int(t) => i for (i, t) in enumerate(node_tags))
    flat_conn = [tag2idx[Int(t)] for t in elem_node_tags[1]]
    tri_nodes = reshape(flat_conn, 3, :)

    gmsh.finalize()
    return xy, tri_nodes
end

"""
    assemble_adjacency(tri_nodes, n_nodes)

Walk every triangle and add a 1 at A[i, j] for each pair (i, j) of its
nodes. The resulting sparse matrix has the same non-zero pattern as a
P1 finite-element stiffness matrix on the same mesh.
"""
function assemble_adjacency(tri_nodes::AbstractMatrix{Int}, n_nodes::Int)
    I = Int[]; J = Int[]; V = Float64[]
    for k in axes(tri_nodes, 2)
        nodes_k = tri_nodes[:, k]              # 3 global node IDs
        for a in 1:3, b in 1:3
            push!(I, nodes_k[a])               # global row
            push!(J, nodes_k[b])               # global column
            push!(V, 1.0)                      # local entry value
        end
    end
    return sparse(I, J, V, n_nodes, n_nodes)   # duplicates auto-summed
end

function fig_mesh_to_matrix()
    @info "Generating Figure 8: Mesh-to-matrix correspondence (real Gmsh.jl run)"

    # Coarse mesh — `lc` chosen so the figure stays readable (~12 nodes).
    xy, tri_nodes = gmsh_unit_square(lc=0.5)
    n_nodes = size(xy, 2)
    n_tri   = size(tri_nodes, 2)

    # Build the matrix exactly as a P1-FEM assembler would.
    A = assemble_adjacency(tri_nodes, n_nodes)

    # Highlight the first triangle so readers can match it to the matrix.
    tri_hl = tri_nodes[:, 1]

    # ---- Left panel: the mesh with node and triangle labels ----
    p1 = plot(legend=false, framestyle=:box, aspect_ratio=:equal,
              xlim=(-0.15, 1.15), ylim=(-0.15, 1.15),
              xlabel="x", ylabel="y",
              title="Mesh — $n_nodes nodes, $n_tri triangles  (T1 highlighted)",
              titlefontsize=11,
              titlefontcolor=parse(Colorant, BLUE))

    for k in 1:n_tri
        a, b, c = tri_nodes[1, k], tri_nodes[2, k], tri_nodes[3, k]
        xs = [xy[1, a], xy[1, b], xy[1, c], xy[1, a]]
        ys = [xy[2, a], xy[2, b], xy[2, c], xy[2, a]]
        fc = (k == 1) ? :gold : :lightblue
        fa = (k == 1) ? 0.55 : 0.20
        plot!(p1, xs, ys, seriestype=:shape,
              fillcolor=fc, fillalpha=fa,
              linecolor=parse(Colorant, BLUE), linewidth=1.4)
        cx = (xy[1, a] + xy[1, b] + xy[1, c]) / 3
        cy = (xy[2, a] + xy[2, b] + xy[2, c]) / 3
        annotate!(p1, cx, cy, text("T$k", :gray30, :center, 8, :italic))
    end

    for i in 1:n_nodes
        scatter!(p1, [xy[1, i]], [xy[2, i]], ms=12, color=:white,
                 markerstrokecolor=parse(Colorant, BLUE), markerstrokewidth=2)
        annotate!(p1, xy[1, i], xy[2, i],
                  text("$i", parse(Colorant, BLUE), :center, 9, :bold))
    end

    # ---- Right panel: spy-style plot of A ----
    rows, cols, _ = findnz(A)

    p2 = scatter(cols, rows, ms=8, color=parse(Colorant, ACCENT),
                 markerstrokecolor=parse(Colorant, BLUE), markerstrokewidth=1,
                 legend=false, framestyle=:box,
                 xlim=(0.5, n_nodes + 0.5), ylim=(0.5, n_nodes + 0.5),
                 yflip=true, aspect_ratio=:equal,
                 xticks=1:n_nodes, yticks=1:n_nodes,
                 xlabel="column j  (node j)", ylabel="row i  (node i)",
                 title="Sparsity of A — red rings = T1's 9 entries",
                 titlefontsize=11,
                 titlefontcolor=parse(Colorant, BLUE))

    # ring the 9 entries contributed by triangle T1
    for i in tri_hl, j in tri_hl
        scatter!(p2, [j], [i], ms=14, color=:transparent,
                 markerstrokecolor=:firebrick, markerstrokewidth=2.5)
    end

    plt = plot(p1, p2, layout=(1, 2),
               size=(CFG.fig_width, 560), dpi=CFG.fig_dpi,
               left_margin=4Plots.mm, bottom_margin=4Plots.mm,
               top_margin=8Plots.mm,
               plot_title="From Gmsh mesh to sparse matrix:  T1 = nodes {$(tri_hl[1]), $(tri_hl[2]), $(tri_hl[3])} ⇒ 9 entries A[i,j] with i,j ∈ T1",
               plot_titlefontsize=11,
               plot_titlefontcolor=parse(Colorant, BLUE))
    savefig(plt, joinpath(CFG.figures_dir, "mesh_to_matrix.png"))
end

"""
    gmsh_geometry(geom_name::Symbol)

Mesh one of four classic 2-D geometries with Gmsh and return
`(xy, tri_nodes)` in the same dense (1..n_nodes) indexing convention as
`gmsh_unit_square`. Each branch demonstrates a different Gmsh feature:

- `:square_refined` — built-in geo kernel with per-point size hints
- `:lshape`         — built-in geo kernel, six-corner re-entrant domain
- `:disc`           — OpenCASCADE kernel, single-call disc primitive
- `:square_hole`    — OpenCASCADE boolean cut (rectangle minus disc)
"""
function gmsh_geometry(geom_name::Symbol)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(string(geom_name))

    if geom_name == :square_refined
        # Same square as before, but the top-left corner asks for a
        # MUCH finer local size — Gmsh smoothly grades from coarse
        # (lc=0.4) to fine (lc=0.05) across the domain.
        p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, 0.4)
        p2 = gmsh.model.geo.addPoint(1.0, 0.0, 0.0, 0.4)
        p3 = gmsh.model.geo.addPoint(1.0, 1.0, 0.0, 0.4)
        p4 = gmsh.model.geo.addPoint(0.0, 1.0, 0.0, 0.05)   # fine corner
        l1 = gmsh.model.geo.addLine(p1, p2)
        l2 = gmsh.model.geo.addLine(p2, p3)
        l3 = gmsh.model.geo.addLine(p3, p4)
        l4 = gmsh.model.geo.addLine(p4, p1)
        cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
        gmsh.model.geo.addPlaneSurface([cl])
        gmsh.model.geo.synchronize()

    elseif geom_name == :lshape
        # Classic re-entrant-corner domain — the hard case for elliptic
        # PDEs because the solution has a singular gradient at the
        # inner corner. Six points walked counter-clockwise.
        coords = [(0.0, 0.0), (1.0, 0.0), (1.0, 0.5),
                  (0.5, 0.5), (0.5, 1.0), (0.0, 1.0)]
        pts = [gmsh.model.geo.addPoint(x, y, 0.0, 0.12) for (x, y) in coords]
        lines = [gmsh.model.geo.addLine(pts[i], pts[mod1(i + 1, length(pts))])
                 for i in eachindex(pts)]
        cl = gmsh.model.geo.addCurveLoop(lines)
        gmsh.model.geo.addPlaneSurface([cl])
        gmsh.model.geo.synchronize()

    elseif geom_name == :disc
        # OpenCASCADE makes curved primitives trivial: one call, done.
        gmsh.model.occ.addDisk(0.0, 0.0, 0.0, 1.0, 1.0)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.18)

    elseif geom_name == :square_hole
        # Rectangle MINUS a smaller disc — boolean operations are why
        # OpenCASCADE exists. (dim, tag) tuples identify entities.
        rect = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 1.0, 1.0)
        hole = gmsh.model.occ.addDisk(0.5, 0.5, 0.0, 0.22, 0.22)
        gmsh.model.occ.cut([(2, rect)], [(2, hole)])
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.12)

    else
        gmsh.finalize()
        error("Unknown geometry: $geom_name")
    end

    gmsh.model.mesh.generate(2)

    node_tags, coord, _   = gmsh.model.mesh.getNodes()
    _, _, elem_node_tags  = gmsh.model.mesh.getElements(2)

    n_nodes = length(node_tags)
    xy      = reshape(coord, 3, n_nodes)[1:2, :]
    tag2idx = Dict(Int(t) => i for (i, t) in enumerate(node_tags))
    flat_conn = [tag2idx[Int(t)] for t in elem_node_tags[1]]
    tri_nodes = reshape(flat_conn, 3, :)

    gmsh.finalize()
    return xy, tri_nodes
end

function fig_gmsh_gallery()
    @info "Generating Figure 9: Gmsh gallery of classic geometries"

    geoms = [:square_refined, :lshape, :disc, :square_hole]
    titles = ["Refined square (fine top-left)", "L-shape (re-entrant corner)",
              "Unit disc (OpenCASCADE)",        "Square with circular hole"]
    xlims = [(-0.10, 1.10), (-0.10, 1.10), (-1.15, 1.15), (-0.10, 1.10)]
    ylims = [(-0.10, 1.10), (-0.10, 1.10), (-1.15, 1.15), (-0.10, 1.10)]

    panels = Plots.Plot[]
    for (i, g) in enumerate(geoms)
        xy, tri_nodes = gmsh_geometry(g)
        n_nodes = size(xy, 2)
        n_tri   = size(tri_nodes, 2)

        p = plot(legend=false, framestyle=:box, aspect_ratio=:equal,
                 xlim=xlims[i], ylim=ylims[i],
                 xlabel="x", ylabel="y",
                 title="$(titles[i]) — $n_nodes nodes, $n_tri tris",
                 titlefontsize=10,
                 titlefontcolor=parse(Colorant, BLUE))

        for k in 1:n_tri
            a, b, c = tri_nodes[1, k], tri_nodes[2, k], tri_nodes[3, k]
            xs = [xy[1, a], xy[1, b], xy[1, c], xy[1, a]]
            ys = [xy[2, a], xy[2, b], xy[2, c], xy[2, a]]
            plot!(p, xs, ys, seriestype=:shape,
                  fillcolor=:lightblue, fillalpha=0.30,
                  linecolor=parse(Colorant, BLUE), linewidth=0.8)
        end

        push!(panels, p)
    end

    plt = plot(panels..., layout=(2, 2),
               size=(CFG.fig_width, 920), dpi=CFG.fig_dpi,
               left_margin=4Plots.mm, bottom_margin=4Plots.mm,
               top_margin=6Plots.mm,
               plot_title="Gmsh recipes — four classic 2-D geometries",
               plot_titlefontsize=13,
               plot_titlefontcolor=parse(Colorant, BLUE))
    savefig(plt, joinpath(CFG.figures_dir, "gmsh_gallery.png"))
end

# =============================================================================
# MAKIE FIGURES  — publication-quality PDE plots rendered headless via CairoMakie
#
# These power the "Plotting PDEs with Makie" section of the website. Every Makie
# call is qualified with `CM.` (see the `import CairoMakie as CM` at the top) so
# the names never clash with the identically-named Plots functions above.
# =============================================================================

# Analytic 3-D scalar field u(x,y,z) = sin(πx) sin(πy) sin(πz) — the exact
# solution of the 3-D Poisson test problem. Returned as (xs, ys, zs, V).
function poisson_field_3d(; N=40, L=CFG.poisson_L)
    h  = L / (N + 1)
    # Keep these as ranges (not collected): Makie's volume/iso-surface path reads
    # the x/y/z arguments as interval endpoints and rejects arbitrary Vectors.
    xs = range(h, L - h, length=N)
    ys = range(h, L - h, length=N)
    zs = range(h, L - h, length=N)
    V  = [sin(π*x) * sin(π*y) * sin(π*z) for x in xs, y in ys, z in zs]
    return xs, ys, zs, V
end

"""
    fig_makie_heatmap_2d()

2-D view of a PDE field: a Makie `heatmap` of the Poisson solution u(x,y) with
an attached colorbar. The clearest first look at any 2-D scalar field.
"""
function fig_makie_heatmap_2d()
    @info "Generating Makie figure: 2-D heatmap"
    xs, ys, U, _ = solve_poisson_2d()

    fig = CM.Figure(size=(760, 620), fontsize=16)
    ax  = CM.Axis(fig[1, 1];
                  title="2-D heatmap — Poisson solution u(x,y)",
                  titlecolor=BLUE, xlabel="x", ylabel="y",
                  aspect=CM.DataAspect())
    hm  = CM.heatmap!(ax, xs, ys, U; colormap=:viridis)
    CM.Colorbar(fig[1, 2], hm; label="u(x,y)")
    CM.save(joinpath(CFG.figures_dir, "makie_heatmap_2d.png"), fig; px_per_unit=2)
end

"""
    fig_makie_surface_3d()

3-D surface of a 2-D PDE field: the same Poisson solution lifted to a Makie
`Axis3` surface, ideal for reading off peaks, valleys and saddles.
"""
function fig_makie_surface_3d()
    @info "Generating Makie figure: 3-D surface"
    xs, ys, U, _ = solve_poisson_2d()

    fig = CM.Figure(size=(820, 640), fontsize=16)
    ax  = CM.Axis3(fig[1, 1];
                   title="3-D surface — u(x,y) over the unit square",
                   titlecolor=BLUE, xlabel="x", ylabel="y", zlabel="u",
                   azimuth=0.6π, elevation=0.18π)
    sp  = CM.surface!(ax, xs, ys, U; colormap=:viridis)
    CM.Colorbar(fig[1, 2], sp; label="u(x,y)")
    CM.save(joinpath(CFG.figures_dir, "makie_surface_3d.png"), fig; px_per_unit=2)
end

"""
    fig_makie_contour()

Contour views side by side: a filled `contourf` (left) and a filled map with
black `contour` iso-lines overlaid (right) — the standard way to show level
sets of a 2-D solution.
"""
function fig_makie_contour()
    @info "Generating Makie figure: filled + line contours"
    xs, ys, U, _ = solve_poisson_2d()
    levels = range(0, maximum(U); length=12)

    fig = CM.Figure(size=(1180, 560), fontsize=16)

    ax1 = CM.Axis(fig[1, 1]; title="Filled contour (contourf)", titlecolor=BLUE,
                  xlabel="x", ylabel="y", aspect=CM.DataAspect())
    cf  = CM.contourf!(ax1, xs, ys, U; levels=levels, colormap=:viridis)
    CM.Colorbar(fig[1, 2], cf; label="u(x,y)")

    ax2 = CM.Axis(fig[1, 3]; title="Filled + iso-lines (contour)", titlecolor=BLUE,
                  xlabel="x", ylabel="y", aspect=CM.DataAspect())
    cf2 = CM.contourf!(ax2, xs, ys, U; levels=levels, colormap=:viridis)
    CM.contour!(ax2, xs, ys, U; levels=levels, color=:black, linewidth=0.9)
    CM.Colorbar(fig[1, 4], cf2; label="u(x,y)")

    CM.save(joinpath(CFG.figures_dir, "makie_contour.png"), fig; px_per_unit=2)
end

"""
    fig_makie_volume_3d()

Visualising a genuinely 3-D field: several z-slices of
u(x,y,z)=sin(πx)sin(πy)sin(πz) drawn as colored `surface` planes stacked inside
an `Axis3`. This is the CairoMakie-friendly way to look into a volume — true
iso-surface/volume rendering needs GLMakie.
"""
function fig_makie_volume_3d()
    @info "Generating Makie figure: 3-D stacked slices"
    xs, ys, zs, V = poisson_field_3d()
    N = length(zs)

    fig = CM.Figure(size=(820, 660), fontsize=16)
    ax  = CM.Axis3(fig[1, 1];
                   title="3-D field — z-slices of u(x,y,z)",
                   titlecolor=BLUE, xlabel="x", ylabel="y", zlabel="z",
                   azimuth=0.55π, elevation=0.16π)

    local sp
    for frac in (0.25, 0.5, 0.75)               # three horizontal cut planes
        k      = round(Int, frac * N)
        zplane = fill(zs[k], length(xs), length(ys))   # flat plane at height z
        sp = CM.surface!(ax, xs, ys, zplane;
                         color=V[:, :, k], colormap=:viridis,
                         colorrange=(0.0, 1.0), transparency=true)
    end
    CM.Colorbar(fig[1, 2], sp; label="u(x,y,z)")
    CM.save(joinpath(CFG.figures_dir, "makie_volume_3d.png"), fig; px_per_unit=2)
end

function generate_all_figures()
    @info "=== Generating all figures ==="
    fig_pde_pipeline()
    fig_pde_types()
    fig_heat_equation()
    fig_wave_equation()
    fig_poisson_equation()
    fig_navier_stokes()
    fig_convergence_study()
    fig_mesh_to_matrix()
    fig_gmsh_gallery()
    # Makie-rendered figures for the "Plotting PDEs with Makie" section
    fig_makie_heatmap_2d()
    fig_makie_surface_3d()
    fig_makie_contour()
    fig_makie_volume_3d()
    @info "All figures saved to $(CFG.figures_dir)"
end

# =============================================================================
# DOCX GENERATION  — uses python-docx via PythonCall
# =============================================================================

"""
Loads python-docx (auto-installs via pip if missing) and returns the module.
"""
function load_python_docx()
    try
        return pyimport("docx")
    catch
        @info "python-docx not found — installing via CondaPkg..."
        # PythonCall's bundled conda env has no pip; use CondaPkg's pip-side helper.
        CondaPkg = Base.require(Base.PkgId(
            Base.UUID("992eb4ea-22a4-4c89-a5bb-47a3300528ab"), "CondaPkg"))
        Base.invokelatest(CondaPkg.add_pip, "python-docx")
        Base.invokelatest(CondaPkg.resolve)
        return pyimport("docx")
    end
end

function build_docx()
    @info "=== Building DOCX ==="
    docx = load_python_docx()
    shared  = pyimport("docx.shared")
    enum_t  = pyimport("docx.enum.text")
    enum_tb = pyimport("docx.enum.table")
    oxml_ns = pyimport("docx.oxml.ns")
    oxml    = pyimport("docx.oxml")

    Inches = shared.Inches
    Pt     = shared.Pt
    RGBColor = shared.RGBColor
    WD_ALIGN = enum_t.WD_ALIGN_PARAGRAPH

    doc = docx.Document()

    # Set default font
    style = doc.styles["Normal"]
    style.font.name = "Arial"
    style.font.size = Pt(11)

    # Page margins
    for section in doc.sections
        section.left_margin   = Inches(1)
        section.right_margin  = Inches(1)
        section.top_margin    = Inches(1)
        section.bottom_margin = Inches(1)
    end

    # ---- helpers -----------------------------------------------------------
    function add_para(text; bold=false, size=11, color=nothing,
                      align=nothing, italic=false)
        p = doc.add_paragraph()
        if align !== nothing; p.alignment = align; end
        run = p.add_run(text)
        run.bold = bold
        run.italic = italic
        run.font.name = "Arial"
        run.font.size = Pt(size)
        if color !== nothing
            run.font.color.rgb = RGBColor(color...)
        end
        return p
    end

    function add_heading(text, level)
        h = doc.add_heading("", level=level)
        run = h.add_run(text)
        run.font.name = "Arial"
        run.bold = true
        if level == 1
            run.font.size = Pt(18)
            run.font.color.rgb = RGBColor(0x1F, 0x4E, 0x79)
        elseif level == 2
            run.font.size = Pt(14)
            run.font.color.rgb = RGBColor(0x2E, 0x75, 0xB6)
        else
            run.font.size = Pt(12)
            run.font.color.rgb = RGBColor(0x44, 0x44, 0x44)
        end
        return h
    end

    function add_code(lines)
        for line in lines
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.25)
            run = p.add_run(line)
            run.font.name = "Courier New"
            run.font.size = Pt(9.5)
            run.font.color.rgb = RGBColor(0x1A, 0x1A, 0x1A)
            # gray background via XML shading
            shading = oxml.OxmlElement("w:shd")
            shading.set(oxml_ns.qn("w:fill"), "F0F0F0")
            p._p.get_or_add_pPr().append(shading)
        end
    end

    function add_bullet(text)
        p = doc.add_paragraph(style="List Bullet")
        run = p.add_run(text)
        run.font.name = "Arial"
        run.font.size = Pt(11)
    end

    function add_image(path; width_in=6.2, caption=nothing)
        if isfile(path)
            # Center via a paragraph that contains the picture
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN.CENTER
            run = p.add_run()
            run.add_picture(path, width=Inches(width_in))
            if caption !== nothing
                add_para(caption; italic=true, size=9, color=(0x55,0x55,0x55),
                         align=WD_ALIGN.CENTER)
            end
        else
            @warn "Image not found: $path"
        end
    end

    function add_table(headers, rows; col_widths=nothing)
        ncol = length(headers)
        tbl = doc.add_table(rows=length(rows)+1, cols=ncol)
        # Try to apply a styled table; fall back to default grid if missing
        try
            tbl.style = "Light Grid Accent 1"
        catch
            try; tbl.style = "Table Grid"; catch; end
        end
        # header row
        for (j, h) in enumerate(headers)
            cell = tbl.rows[0].cells[j-1]
            cell.text = ""
            run = cell.paragraphs[0].add_run(h)
            run.bold = true
            run.font.name = "Arial"
            run.font.size = Pt(10)
            run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
            tcPr = cell._tc.get_or_add_tcPr()
            shd = oxml.OxmlElement("w:shd")
            shd.set(oxml_ns.qn("w:fill"), "1F4E79")
            tcPr.append(shd)
        end
        # data rows
        for (i, row) in enumerate(rows)
            for (j, cell_text) in enumerate(row)
                c = tbl.rows[i].cells[j-1]
                c.text = ""
                r = c.paragraphs[0].add_run(string(cell_text))
                r.font.name = "Arial"
                r.font.size = Pt(10)
            end
        end
        return tbl
    end

    function add_page_break()
        doc.add_page_break()
    end

    function add_toc()
        p_heading = doc.add_paragraph()
        run_h = p_heading.add_run("Table of Contents")
        run_h.bold = true
        run_h.font.name = "Arial"
        run_h.font.size = Pt(22)
        run_h.font.color.rgb = RGBColor(0x1F, 0x4E, 0x79)
        p_heading.alignment = WD_ALIGN.CENTER

        p = doc.add_paragraph()
        run = p.add_run()
        fldChar_begin = oxml.OxmlElement("w:fldChar")
        fldChar_begin.set(oxml_ns.qn("w:fldCharType"), "begin")
        run._r.append(fldChar_begin)
        instrText = oxml.OxmlElement("w:instrText")
        instrText.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        instrText.text = " TOC \\o \"1-3\" \\h \\z \\u "
        run._r.append(instrText)
        fldChar_sep = oxml.OxmlElement("w:fldChar")
        fldChar_sep.set(oxml_ns.qn("w:fldCharType"), "separate")
        run._r.append(fldChar_sep)
        fldChar_end = oxml.OxmlElement("w:fldChar")
        fldChar_end.set(oxml_ns.qn("w:fldCharType"), "end")
        run._r.append(fldChar_end)
    end

    # ─── TITLE PAGE ──────────────────────────────────────────────────────
    for _ in 1:3; doc.add_paragraph(); end
    add_para(CFG.doc_title; bold=true, size=36, color=(0x1F,0x4E,0x79),
             align=WD_ALIGN.CENTER)
    add_para(CFG.doc_subtitle; bold=true, size=22, color=(0x2E,0x75,0xB6),
             align=WD_ALIGN.CENTER)
    add_para(CFG.doc_tagline; size=12, color=(0x66,0x66,0x66),
             align=WD_ALIGN.CENTER)
    for _ in 1:2; doc.add_paragraph(); end
    add_para("Author: $(CFG.doc_author)"; size=11, color=(0x55,0x55,0x55),
             align=WD_ALIGN.CENTER)
    add_para("Generated: $(string(Dates.today()))"; size=10, italic=true,
             color=(0x88,0x88,0x88), align=WD_ALIGN.CENTER)
    add_page_break()

    # ─── TABLE OF CONTENTS ───────────────────────────────────────────────
    add_toc()
    add_page_break()

    # ─── OVERVIEW ───────────────────────────────────────────────────────
    add_heading("Overview & Prerequisites", 1)
    add_para("This guide walks you through setting up a complete Julia " *
             "environment for scientific computing and HPC on Windows. " *
             "Steps 1–3 install Julia. Step 4 introduces PDE packages. " *
             "Step 4.5 builds the juliaPDEs package — the same one shipped " *
             "in the juliaPDEs/ folder of this repo — so you have a real " *
             "module with Revise live-reload and a test suite before " *
             "tackling the worked examples. Step 5 walks the full PDE-" *
             "solving workflow with worked examples for the heat, wave, " *
             "Poisson, and Navier-Stokes equations and a section on " *
             "unstructured meshing with Gmsh. Steps 6–8 connect to the " *
             "HPC cluster.")
    add_page_break()

    # ─── STEP 1: WSL ────────────────────────────────────────────────────
    add_heading("Step 1 — Install WSL2", 1)
    add_para("WSL2 runs a real Linux kernel inside Windows. Open " *
             "PowerShell as Administrator and run:")
    add_code(["wsl --install",
              "# Then restart your computer when prompted"])
    add_para("After restart, Ubuntu opens automatically; create a Linux " *
             "username and password, then update the system:")
    add_code(["sudo apt update && sudo apt upgrade -y",
              "sudo apt install -y curl wget git build-essential"])
    add_page_break()

    # ─── STEP 2: JULIA ──────────────────────────────────────────────────
    add_heading("Step 2 — Install juliaup, then Julia", 1)
    add_para("juliaup is the official Julia version manager — like pyenv " *
             "for Python or nvm for Node.")
    add_heading("2.1  Install juliaup", 2)
    add_code(["curl -fsSL https://install.julialang.org | sh",
              "source ~/.bashrc",
              "juliaup --version"])
    add_heading("2.2  Use juliaup to install Julia", 2)
    add_code(["juliaup add release       # latest stable",
              "juliaup default release   # set as default",
              "julia --version"])
    add_table(["Command", "Action"],
              [["juliaup add release",   "Install latest stable Julia"],
               ["juliaup add 1.10.4",    "Install a specific version"],
               ["juliaup default 1.10.4","Switch the default version"],
               ["juliaup list",          "List installed versions"],
               ["juliaup update",        "Update everything"]])
    add_page_break()

    # ─── STEP 7: VS CODE ────────────────────────────────────────────────
    add_heading("Step 3 — Configure VS Code for Julia + WSL", 1)
    add_para("Install these VS Code extensions:")
    add_table(["Extension", "Purpose"],
              [["WSL",            "Connect VS Code to WSL"],
               ["Julia",          "Syntax, IntelliSense, REPL"],
               ["Jupyter",        "Run notebooks with Julia kernel"],
               ["Remote — SSH",   "Connect to HPC cluster"]])
    add_page_break()

    # ─── STEP 3: PACKAGES ───────────────────────────────────────────────
    add_heading("Step 4 — Installing Julia Packages", 1)
    add_para("Julia's package manager (Pkg) is built into the language. " *
             "Press ] in the REPL to enter package mode.")
    add_heading("3.1  Project Environments (best practice)", 2)
    add_code(["mkdir -p ~/julia-pde && cd ~/julia-pde",
              "julia --project=.",
              "# In package mode (press ]):",
              "(@julia-pde) pkg> add DifferentialEquations Plots SparseArrays"])
    add_para("Julia creates Project.toml and Manifest.toml — commit both to " *
             "git for full reproducibility on the cluster.")
    add_page_break()

    # ─── STEP 4: PDE PACKAGES ───────────────────────────────────────────
    add_heading("Step 5 — Key Julia Packages for PDE Solving", 1)
    add_table(["Package", "Use case"],
              [["LinearAlgebra (stdlib)",     "Dense matrix ops"],
               ["SparseArrays (stdlib)",      "Sparse PDE matrices"],
               ["FFTW.jl",                    "Spectral methods"],
               ["DifferentialEquations.jl",   "ODE/PDE solver suite"],
               ["MethodOfLines.jl",           "Auto-discretization of PDEs"],
               ["KrylovKit.jl",               "GMRES, CG iterative solvers"],
               ["Gmsh.jl",                    "Unstructured mesh generation"],
               ["Ferrite.jl / Gridap.jl",     "Finite element methods"],
               ["Oceananigans.jl",            "Production CFD on HPC"],
               ["Plots.jl / Makie.jl",        "Visualization"],
               ["MPI.jl / CUDA.jl",           "HPC parallelism"]])
    add_page_break()

    # ─── STEP 4.5: BUILD JULIAPDES PACKAGE ──────────────────────────────
    add_heading("Step 6 — Build the juliaPDEs Package (was Step 4.5)", 1)
    add_para("Before tackling the worked examples in Step 5, package the " *
             "code into a real Julia package called juliaPDEs. The same " *
             "module — with src/heat.jl, src/wave.jl, and " *
             "src/Navier_Stokes.jl — is what ships in the juliaPDEs/ folder " *
             "of this repo. Working from a package gives you Revise live-" *
             "reload, a test suite, and one stable import you can call from " *
             "the REPL, from tests, and from Slurm job scripts.")

    add_heading("4.5.1  Step A — Generate the skeleton", 2)
    add_para("One command creates the folder and the boilerplate Project.toml " *
             "plus module file. Add a test/ folder yourself — Julia uses it " *
             "automatically:")
    add_code(["using Pkg",
              "Pkg.generate(\"juliaPDEs\")   # creates juliaPDEs/Project.toml + src/juliaPDEs.jl",
              "mkdir juliaPDEs/test          # add a test/ folder yourself",
              "touch juliaPDEs/test/runtests.jl"])

    add_heading("4.5.2  Step B — Add the dependencies", 2)
    add_para("Activate the new package's own environment and install everything " *
             "the three solvers need. Plots is for visualization, Revise for " *
             "live reload, and Oceananigans is the CFD engine called by the " *
             "Navier-Stokes solver. The remaining packages match the " *
             "Project.toml that ships with the repo:")
    add_code(["cd juliaPDEs",
              "julia --project=.",
              "# In package mode (press ]):",
              "(@juliaPDEs) pkg> add Plots Revise Oceananigans LinearAlgebra SparseArrays Gmsh Makie CUDA"])
    add_para("Julia writes the deps into Project.toml and pins exact versions " *
             "in Manifest.toml — commit both to git so the cluster gets the " *
             "same versions you developed against.")

    add_heading("4.5.3  Step C — The module entry point: src/juliaPDEs.jl", 2)
    add_para("Everything inside `module … end` is private by default. " *
             "`using` brings third-party packages into the module once. " *
             "`include()` pastes a sub-file in so it sees the module's scope. " *
             "`export` declares the names visible after `using juliaPDEs`:")
    add_code(["module juliaPDEs",
              "",
              "using Plots, Revise, Oceananigans",
              "",
              "include(\"heat.jl\")",
              "include(\"wave.jl\")",
              "include(\"Navier_Stokes.jl\")",
              "",
              "export solve_heat_1d, solve_wave_1d, solve_navier_stokes,",
              "       animate_wave_1d, animate_navier_stokes",
              "",
              "end # module juliaPDEs"])
    add_page_break()

    add_heading("4.5.4  Step D — The heat solver: src/heat.jl", 2)
    add_para("A 1-D explicit Euler heat solver: Dirichlet BCs (the endpoints " *
             "stay 0), Gaussian initial pulse. The CFL factor 0.4 keeps " *
             "r = α·dt/dx² ≤ 0.5 so the explicit scheme is stable.")
    add_code(["function solve_heat_1d(; N=200, α=0.01, L=1.0, T=1.0)",
              "    dx     = L / (N + 1)",
              "    x      = range(dx, L - dx, length=N)",
              "    dt     = 0.4 * dx^2 / α          # CFL: r ≤ 0.5",
              "    nsteps = ceil(Int, T / dt)",
              "    u = @. exp(-100 * (x - 0.5)^2)   # Gaussian initial condition",
              "    for _ in 1:nsteps",
              "        # Explicit Euler + 3-point stencil:",
              "        #   u_i^{n+1} = u_i^n + r·(u_{i+1} - 2u_i + u_{i-1})",
              "        u[2:end-1] .+= (α * dt / dx^2) .*",
              "            (u[3:end] .- 2 .* u[2:end-1] .+ u[1:end-2])",
              "        # Dirichlet BCs: u[1] = u[end] = 0 — never updated",
              "    end",
              "    return collect(x), u",
              "end"])
    add_page_break()

    add_heading("Beginner Concept Check: Structs, Instances, and Functions", 3)
    add_para("Notice how the Heat solver above uses plain keyword arguments. As " *
             "simulations grow more complex, passing individual variables becomes " *
             "unmanageable. Before moving to the Wave solver, let's clarify how " *
             "custom data types solve this:")
    add_para("1. What is a Struct? A struct is a custom composite data type that " *
             "bundles multiple pieces of related data (called fields) together under " *
             "a single name. Instead of passing separate variables into every " *
             "function, we pack them into one structured container blueprint.")
    add_para("2. How to create a Struct? We write the struct keyword followed by " *
             "the container name, list its internal fields with their types, and end " *
             "with end. Using the Base.@kwdef macro above it automatically writes " *
             "a keyword constructor with convenient default values:")
    add_code(["Base.@kwdef struct SimpleHeatConfig",
              "    L::Float64 = 1.0       # Domain length",
              "    α::Float64 = 0.01      # Thermal diffusivity",
              "    N::Int     = 100       # Grid resolution",
              "end"])
    add_para("3. What is the variable that is the Struct (an Instance)? The code " *
             "above is just a blueprint. To use it, we call the struct's name like " *
             "a function to create a concrete object in computer memory, called an " *
             "instance, and assign it to a variable. We can then read its internal " *
             "fields using dot notation:")
    add_code(["# Create an instance; omitted fields automatically use defaults",
              "my_config = SimpleHeatConfig(α=0.05, N=200)",
              "",
              "# Access individual bundled fields using a dot (.):",
              "println(\"Configured diffusivity: \", my_config.α)",
              "println(\"Configured grid points: \", my_config.N)"])
    add_para("4. How to create a Function to receive the Struct? To write a function " *
             "that operates on our configuration, we define a function that accepts " *
             "a single argument annotated with our custom struct type. Inside the " *
             "function body, we use dot notation to unpack the exact fields we need:")
    add_code(["function print_stability_limit(config::SimpleHeatConfig)",
              "    dx = config.L / (config.N + 1)",
              "    max_dt = 0.5 * dx^2 / config.α",
              "    println(\"Maximum stable time step dt: \", max_dt)",
              "end",
              "",
              "# Pass our struct variable into the function",
              "print_stability_limit(my_config)"])
    add_para("5. Structs vs. Functions: Nouns vs. Verbs. A struct holds passive state " *
             "and data (e.g., domain dimensions, arrays, physical constants). A " *
             "function defines active behavior and logic (e.g., loops, time-stepping, " *
             "mathematical updates). Structs are the materials; functions are the " *
             "machines that transform them.")
    add_para("6. Global Scope vs. Local Scope (Crucial for Performance). Variables " *
             "defined directly in the main script body live in the global scope. Because " *
             "global variables can change their data type at any time, the Julia " *
             "compiler cannot optimize them, making loops in global scope very slow. " *
             "Variables defined inside a function body live in a local scope. The compiler " *
             "knows their exact types, compiling local operations to extremely fast " *
             "machine code. Best Practice: Always wrap your solver loops and computations " *
             "inside functions!")
    add_para("7. Where do you run these scripts? Once you write your code in a file " *
             "(e.g., solver.jl), you can run it in three ways: (a) Terminal execution: " *
             "run 'julia solver.jl' directly from your Linux command line. (b) REPL " *
             "inclusion: open the Julia REPL and type 'include(\"solver.jl\")' to load " *
             "and execute the file. (c) VS Code execution: open the file in VS Code and " *
             "press Shift+Enter to run specific blocks interactively.")

    add_heading("4.5.5  Step E — The wave solver: src/wave.jl", 2)
    add_para("The wave solver groups its parameters in a WaveEquation struct " *
             "so callers can mix-and-match grid size, boundary condition " *
             "(:Dirichlet, :Neumann, :Periodic, :Absorbing), and initial " *
             "shape independently. Base.@kwdef auto-generates a keyword " *
             "constructor with defaults, so WaveEquation(nx=400, T=2.0) " *
             "leaves everything else default:")
    add_code(["Base.@kwdef struct WaveEquation",
              "    nx::Int = 300",
              "    nt::Int = 300",
              "    c::Float64 = 1.0",
              "    L::Float64 = 1.0",
              "    T::Float64 = 1.5",
              "    bc::Symbol = :Dirichlet  # :Dirichlet | :Neumann | :Periodic | :Absorbing",
              "    f_init::Function = x -> sin.(pi * x)",
              "    dx::Float64 = L / (nx - 1)",
              "    dt::Float64 = T / nt",
              "    x::Any = range(0, L, length=nx)",
              "end"])
    add_para("The leap-frog update uses the same 3-point Laplacian stencil " *
             "as the heat equation, but stepped in time with two history " *
             "slots (u_prev, u_curr). CFL: λ = c·dt/dx ≤ 1.")
    add_code(["function solve_wave_1d(problem::WaveEquation; return_history=false)",
              "    dx, dt, x = problem.dx, problem.dt, problem.x",
              "    λ = (problem.c * dt / dx)^2",
              "    u_curr = problem.f_init.(x)",
              "    apply_boundary_conditions!(u_curr, u_curr, problem; is_initial=true)",
              "    u_prev = copy(u_curr); u_next = copy(u_curr)",
              "    for i in 1:problem.nt",
              "        # Leap-frog: u^{n+1} = 2u^n - u^{n-1} + λ²·(stencil)",
              "        u_next[2:end-1] .= 2 .* u_curr[2:end-1] .- u_prev[2:end-1] .+",
              "            λ .* (u_curr[3:end] .- 2 .* u_curr[2:end-1] .+ u_curr[1:end-2])",
              "        apply_boundary_conditions!(u_next, u_curr, problem)",
              "        u_prev, u_curr = copy(u_curr), copy(u_next)",
              "    end",
              "    return collect(x), u_curr",
              "end"])
    add_para("The full file also defines apply_boundary_conditions! (one " *
             "branch per BC symbol), animate_wave_1d (saves a GIF), and an " *
             "analytic wave_1d_exact reference solution — see " *
             "juliaPDEs/src/wave.jl in the repo for the complete code.")
    add_page_break()

    add_heading("4.5.6  Step F — The Navier-Stokes solver: src/Navier_Stokes.jl", 2)
    add_para("The incompressible Navier-Stokes equations are two coupled " *
             "PDEs in the velocity field u(x,y,t) = (u, v) and the pressure " *
             "p(x,y,t):")
    add_code(["∂u/∂t + (u · ∇)u = -∇p + ν·∇²u        # momentum",
              "             ∇ · u = 0                # incompressibility"])
    add_para("Each term has a physical meaning: (u·∇)u is advection (fluid " *
             "carrying its own momentum — the nonlinear part), -∇p is the " *
             "pressure gradient that enforces incompressibility, ν·∇²u is " *
             "viscous diffusion, and ∇·u = 0 keeps the flow divergence-free " *
             "at every step.")
    add_para("We solve the classic lid-driven cavity benchmark on Ω = " *
             "[0,1]×[0,1]:")
    add_bullet("Top lid (north):  u = 1, v = 0 — drives the flow.")
    add_bullet("Other walls:      u = v = 0 (no-slip) — the default.")
    add_bullet("Initial condition: u(x,y,0) = 0.")
    add_bullet("Reynolds number:  Re = U·L / ν = (1·1)/0.01 = 100, set by " *
               "choosing ν = 1e-2 in ScalarDiffusivity.")
    add_para("Rather than discretising and writing a projection method by " *
             "hand, we delegate to Oceananigans.jl. Oceananigans handles the " *
             "staggered finite-volume grid, a 5th-order upwind-biased " *
             "advection scheme, the divergence-free pressure projection, and " *
             "a 3rd-order Runge-Kutta time integrator — all on CPU or GPU. " *
             "The simulation runs for t ∈ [0, 10] with Δt = 1e-3 and pushes " *
             "a snapshot of the speed field √(u² + v²) into speed_hist every " *
             "100 iterations so we can animate it later with " *
             "animate_navier_stokes:")
    add_code(["function solve_navier_stokes()",
              "    grid = Oceananigans.RectilinearGrid(size=(64, 64),",
              "        x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))",
              "",
              "    # Top lid moves at speed 1",
              "    u_bcs = Oceananigans.FieldBoundaryConditions(",
              "        north = Oceananigans.ValueBoundaryCondition(1.0))",
              "",
              "    model = Oceananigans.NonhydrostaticModel(grid;",
              "        advection   = Oceananigans.UpwindBiased(order=5),",
              "        timestepper = :RungeKutta3,",
              "        boundary_conditions = (u=u_bcs,),",
              "        closure     = Oceananigans.ScalarDiffusivity(ν=1e-2))   # ν = 1/Re",
              "",
              "    Oceananigans.set!(model, u=0, v=0)",
              "    simulation = Oceananigans.Simulation(model, Δt=1e-3, stop_time=10.0)",
              "",
              "    speed_field = Oceananigans.Field(",
              "        sqrt(model.velocities.u^2 + model.velocities.v^2))",
              "    speed_hist  = []",
              "",
              "    function save_state(sim)",
              "        Oceananigans.compute!(speed_field)",
              "        push!(speed_hist, copy(Oceananigans.interior(speed_field, :, :, 1)))",
              "    end",
              "    simulation.callbacks[:save] = Oceananigans.Callback(",
              "        save_state, Oceananigans.IterationInterval(100))",
              "",
              "    Oceananigans.run!(simulation)",
              "    return (speed_hist=speed_hist, Δt=0.1)",
              "end"])
    add_page_break()

    add_heading("4.5.7  Step G — Load the package & run it", 2)
    add_para("Pkg.develop links the package into your environment by path — " *
             "no install, no version pin. Edits to src/ files take effect " *
             "immediately (with Revise, see 4.5.9):")
    add_code(["# From the directory that CONTAINS the juliaPDEs/ folder:",
              "using Pkg",
              "Pkg.develop(path=\"./juliaPDEs\")",
              "",
              "using juliaPDEs",
              "x, u  = solve_heat_1d()                            # heat",
              "x, u  = solve_wave_1d(juliaPDEs.WaveEquation())    # wave (defaults)",
              "hist  = solve_navier_stokes()                      # Navier-Stokes",
              "animate_navier_stokes(hist)                        # → .gif"])

    add_heading("4.5.8  Step H — Write tests", 2)
    add_para("Tests live in test/runtests.jl. Run them with `(@juliaPDEs) " *
             "pkg> test`:")
    add_code(["using juliaPDEs, Test",
              "",
              "@testset \"Heat equation\" begin",
              "    x, u = solve_heat_1d(N=50, T=0.1)",
              "    @test length(x) == 50",
              "    @test maximum(u) < 1.0",
              "    @test u[1]   ≈ 0.0 atol=1e-10",
              "    @test u[end] ≈ 0.0 atol=1e-10",
              "end",
              "",
              "@testset \"Wave equation\" begin",
              "    prob = juliaPDEs.WaveEquation(nx=100, nt=200, T=0.5)",
              "    x, u = solve_wave_1d(prob)",
              "    @test length(x) == 100",
              "    @test all(isfinite, u)",
              "end"])

    add_heading("4.5.9  Revise — live reload during development", 2)
    add_para("Add Revise to your global environment so changes to src/ files " *
             "take effect without restarting Julia:")
    add_code(["(@v1.10) pkg> add Revise"])
    add_para("Make Revise auto-load in every session by adding this to " *
             "~/.julia/config/startup.jl:")
    add_code(["try",
              "    using Revise",
              "catch e",
              "    @warn \"Could not load Revise\" exception = e",
              "end"])
    add_para("Then your day-to-day loop is: edit a src/*.jl file, save, call " *
             "the function in the REPL — the new version runs immediately. " *
             "Revise handles function bodies and module-level code; struct " *
             "redefinitions and changes to top-level `using` still need a " *
             "restart.")
    add_page_break()

    # ─── STEP 5: PDE WORKFLOW ───────────────────────────────────────────
    add_heading("Step 7 — The PDE Solving Workflow (was Step 5)", 1)
    add_para("Solving any PDE follows the same five-stage pipeline. " *
             "Errors at early stages cascade through everything that follows.")

    add_image(joinpath(CFG.figures_dir, "pde_pipeline.png"); width_in=6.5,
              caption="Figure 1. The five-stage PDE-solving pipeline.")
    add_page_break()

    add_heading("5.1  Stage 1 — The Mathematical Model", 2)
    add_para("Every complete PDE problem has four parts: the equation, " *
             "the domain, boundary conditions, and (for time-dependent " *
             "problems) initial conditions.")
    add_image(joinpath(CFG.figures_dir, "pde_types.png"); width_in=6.5,
              caption="Figure 2. Four major PDE families.")
    add_page_break()

    add_heading("5.2  Stage 2 — Discretization", 2)
    add_para("Replace continuous derivatives with finite-difference stencils:")
    add_code(["# 2nd derivative (central, 2nd-order):",
              "#   d²u/dx² ≈ ( u[i+1] - 2u[i] + u[i-1] ) / h²",
              "# 2D Laplacian (5-point stencil):",
              "#   ∇²u ≈ ( u[i+1,j]+u[i-1,j]+u[i,j+1]+u[i,j-1]-4u[i,j] ) / h²"])
    add_para("CFL stability for explicit schemes:")
    add_code(["# Heat (1D, explicit Euler): dt ≤ dx² / (2α)",
              "# Wave (1D):                  dt ≤ dx / c"])

    add_heading("5.3  Stage 3 — Build the System", 2)
    add_code(["using SparseArrays, LinearAlgebra",
              "function laplacian_2d(N, h)",
              "    e   = ones(N)",
              "    L1D = spdiagm(0 => -2*e, 1 => e[1:end-1], -1 => e[1:end-1])",
              "    I_N = sparse(I, N, N)",
              "    return (kron(I_N, L1D) + kron(L1D, I_N)) / h^2",
              "end"])
    add_page_break()

    add_heading("5.4  Stage 4 — Choose & Run a Solver", 2)
    add_table(["Situation", "Solver"],
              [["Elliptic, N < 10⁴",          "u = A \\ b (sparse direct)"],
               ["Elliptic, large symmetric",  "CG (IterativeSolvers.jl)"],
               ["Elliptic, large non-sym.",   "GMRES + ILU"],
               ["Parabolic, small N",         "Explicit Euler"],
               ["Parabolic, stiff",           "Crank-Nicolson / Rosenbrock23"],
               ["Hyperbolic",                 "Leap-frog / Tsit5"],
               ["Navier-Stokes",              "Oceananigans.jl / projection"]])

    # ── Worked Example A: Heat ──
    add_heading("Example A — Heat Equation (parabolic)", 3)
    add_para(@sprintf("Parameters: N=%d, α=%.3f, L=%.1f, T=%.1f",
                      CFG.heat_N, CFG.heat_alpha, CFG.heat_L, CFG.heat_T_final))
    add_image(joinpath(CFG.figures_dir, "heat_equation.png"); width_in=6.5,
              caption="Figure 3. Heat equation: pulse diffuses and smooths over time.")
    add_page_break()

    # ── Worked Example B: Wave ──
    add_heading("Example B — Wave Equation (hyperbolic)", 3)
    add_para(@sprintf("Parameters: N=%d, c=%.1f, L=%.1f, T=%.1f",
                      CFG.wave_N, CFG.wave_c, CFG.wave_L, CFG.wave_T_final))
    add_image(joinpath(CFG.figures_dir, "wave_equation.png"); width_in=6.5,
              caption="Figure 4. Wave equation: pluck splits and reflects off boundaries.")
    add_page_break()

    # ── Worked Example C: Poisson ──
    add_heading("Example C — Poisson Equation (elliptic)", 3)
    add_para(@sprintf("Parameters: N=%d (interior), L=%.1f, single sparse direct solve",
                      CFG.poisson_N, CFG.poisson_L))
    add_image(joinpath(CFG.figures_dir, "poisson_equation.png"); width_in=6.5,
              caption="Figure 5. Poisson equation: source f → solution u, no time evolution.")
    add_page_break()

    # ── Worked Example D: Navier-Stokes ──
    add_heading("Example D — Navier-Stokes (the hard one)", 3)
    add_para("The Navier-Stokes equations describe how fluids flow. They " *
             "are nonlinear, couple velocity and pressure, and admit " *
             "turbulent solutions at high Reynolds numbers.")
    add_para("Equations (incompressible):")
    add_code(["  ∂u/∂t + (u·∇)u  =  -∇p/ρ + ν∇²u + g     ← momentum",
              "       ∇·u           =  0                       ← incompressibility"])
    add_table(["Term", "Meaning"],
              [["∂u/∂t",      "Local acceleration"],
               ["(u·∇)u",     "Convective acceleration (nonlinear!)"],
               ["-∇p/ρ",      "Pressure gradient force"],
               ["ν∇²u",       "Viscous diffusion (friction)"],
               ["∇·u = 0",    "Mass conservation / incompressibility"]])
    add_image(joinpath(CFG.figures_dir, "navier_stokes.png"); width_in=6.5,
              caption=@sprintf("Figure 6. Lid-driven cavity at Re=%.0f. Lid drives a primary vortex.",
                               CFG.ns_Re))
    add_para("Why it's hard:")
    add_bullet("Nonlinear convective term (u·∇)u — no linear superposition.")
    add_bullet("Velocity-pressure coupling — pressure has no time derivative; " *
               "must enforce ∇·u = 0 via projection methods.")
    add_bullet("Turbulence at high Re requires very fine grids or RANS/LES/DNS.")
    add_bullet("Millennium Prize problem — smoothness in 3D is unproven!")
    add_para("Recommended: use Oceananigans.jl for HPC fluid simulations — " *
             "it is Slurm + MPI + GPU ready out of the box.")
    add_page_break()

    add_heading("5.5  Stage 5 — Visualize & Verify", 2)
    add_para("Three levels of checks every PDE solver should pass:")
    add_bullet("Sanity: plot the solution; check BCs are satisfied.")
    add_bullet("Convergence study: refine grid, check error scales as O(h²).")
    add_bullet("Physical validation: compare against analytical solution or benchmark.")
    add_image(joinpath(CFG.figures_dir, "convergence_study.png"); width_in=6.5,
              caption="Figure 7. Left: Poisson L∞ error vs grid spacing h on a log-log " *
                      "scale — the slope confirms 2nd-order convergence. Right: solve " *
                      "time scales near-linearly in N² (sparse direct solver).")

    add_page_break()

    # ─── 5.6  GMSH: MESH → MATRIX ───────────────────────────────────────
    add_heading("5.6  Beyond Structured Grids — Meshing with Gmsh", 2)
    add_para("Finite-difference grids work for rectangles. Real geometry " *
             "— L-shapes, airfoils, brain MRI, machine parts — needs " *
             "unstructured meshes. Gmsh is the de-facto open-source mesher; " *
             "the Julia bindings live in Gmsh.jl. The workflow has three " *
             "parts: (1) describe the geometry, (2) let Gmsh triangulate it, " *
             "(3) read the node coordinates and element connectivity back " *
             "into Julia and assemble a sparse matrix.")

    add_heading("Install Gmsh.jl", 3)
    add_code(["# In package mode (press ]):",
              "(@julia-pde) pkg> add Gmsh"])

    add_heading("Generate a 2D triangular mesh of the unit square", 3)
    add_para("This is the exact code that produced Figure 8 — running the " *
             "lecture script reproduces the same mesh on the student's " *
             "machine, byte-for-byte.")
    add_code(["using Gmsh: gmsh",
              "using SparseArrays",
              "",
              "function gmsh_unit_square(; lc::Float64 = 0.5)",
              "    gmsh.initialize()",
              "    gmsh.option.setNumber(\"General.Terminal\", 0)",
              "    gmsh.model.add(\"unit_square\")",
              "",
              "    # 1. Geometry — corner points of [0,1]×[0,1], target size lc",
              "    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)",
              "    p2 = gmsh.model.geo.addPoint(1.0, 0.0, 0.0, lc)",
              "    p3 = gmsh.model.geo.addPoint(1.0, 1.0, 0.0, lc)",
              "    p4 = gmsh.model.geo.addPoint(0.0, 1.0, 0.0, lc)",
              "",
              "    # 2. Edges and the surface they bound",
              "    l1 = gmsh.model.geo.addLine(p1, p2)",
              "    l2 = gmsh.model.geo.addLine(p2, p3)",
              "    l3 = gmsh.model.geo.addLine(p3, p4)",
              "    l4 = gmsh.model.geo.addLine(p4, p1)",
              "    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])",
              "    gmsh.model.geo.addPlaneSurface([cl])",
              "",
              "    # 3. Synchronise the CAD kernel and triangulate the surface",
              "    gmsh.model.geo.synchronize()",
              "    gmsh.model.mesh.generate(2)        # 2D ⇒ triangles",
              "",
              "    # 4. Pull the mesh back into Julia",
              "    node_tags, coord, _  = gmsh.model.mesh.getNodes()",
              "    _, _, elem_node_tags = gmsh.model.mesh.getElements(2)",
              "",
              "    n_nodes = length(node_tags)",
              "    xy      = reshape(coord, 3, n_nodes)[1:2, :]   # drop z",
              "",
              "    # Gmsh tags can be sparse — remap them to dense 1..n_nodes",
              "    tag2idx   = Dict(Int(t) => i for (i,t) in enumerate(node_tags))",
              "    flat_conn = [tag2idx[Int(t)] for t in elem_node_tags[1]]",
              "    tri_nodes = reshape(flat_conn, 3, :)",
              "",
              "    gmsh.finalize()",
              "    return xy, tri_nodes",
              "end",
              "",
              "xy, tri_nodes = gmsh_unit_square(lc=0.5)",
              "# xy        :: 2 × n_nodes  (coordinates of every mesh node)",
              "# tri_nodes :: 3 × n_tri    (column k = global node IDs of triangle k)"])

    add_heading("From mesh to sparse matrix — element assembly", 3)
    add_para("Visit every triangle, and for each pair (a, b) of its three " *
             "local nodes add the local contribution at the matching global " *
             "indices. With unit local entries this gives a matrix whose " *
             "non-zero pattern is identical to a P1 stiffness matrix on the " *
             "same mesh:")
    add_code(["function assemble_adjacency(tri_nodes, n_nodes)",
              "    I = Int[]; J = Int[]; V = Float64[]",
              "    for k in axes(tri_nodes, 2)",
              "        nodes_k = tri_nodes[:, k]      # 3 global node IDs",
              "        for a in 1:3, b in 1:3",
              "            push!(I, nodes_k[a])       # global row",
              "            push!(J, nodes_k[b])       # global column",
              "            push!(V, 1.0)              # local entry value",
              "        end",
              "    end",
              "    return sparse(I, J, V, n_nodes, n_nodes)  # duplicates auto-summed",
              "end",
              "",
              "n_nodes = size(xy, 2)",
              "A = assemble_adjacency(tri_nodes, n_nodes)"])
    add_para("The COO triplets (I, J, V) are the natural output of element " *
             "assembly; `sparse` automatically sums duplicate (i, j) pairs, " *
             "which is exactly what \"adding contributions from every " *
             "element touching that node\" means. Replace the `1.0` with " *
             "the entries of a 3×3 local stiffness matrix to get a real FEM " *
             "Laplacian — the index plumbing is the same.")

    add_heading("How mesh indices and matrix indices line up", 3)
    add_para("This is the central idea behind FEM/FVM assembly — once it " *
             "clicks, every other implementation detail follows:")
    add_bullet("Each mesh node has a unique global index 1..n. " *
               "That same index is its row (and column) in the matrix A.")
    add_bullet("Each element carries a small connectivity list. For a P1 " *
               "triangle that list has three node IDs.")
    add_bullet("The element produces a small dense local matrix whose rows " *
               "and columns are labelled by local node order (1, 2, 3).")
    add_bullet("Local index a maps to global index nodes_k[a], so " *
               "K_local[a, b] is added into A[nodes_k[a], nodes_k[b]].")
    add_bullet("If two elements share a node, both contribute to the same " *
               "row/column of A — A[i, j] is the sum of every element " *
               "contribution touching the pair (i, j). Mesh edges become " *
               "off-diagonal non-zeros; mesh nodes become diagonal entries.")

    add_image(joinpath(CFG.figures_dir, "mesh_to_matrix.png"); width_in=6.5,
              caption="Figure 8. Left: a 9-node, 8-triangle mesh of the unit " *
                      "square. Right: the sparsity pattern of the assembled " *
                      "matrix A. An entry A[i, j] is non-zero iff nodes i " *
                      "and j share at least one triangle. The red rings on " *
                      "the right mark the nine entries contributed by " *
                      "triangle T1 (gold on the left), whose nodes are " *
                      "{1, 2, 5} — those nine entries are A[1,1], A[1,2], " *
                      "A[1,5], A[2,1], A[2,2], A[2,5], A[5,1], A[5,2], A[5,5].")

    add_para("Practical implications:")
    add_bullet("Sparsity ≈ average node degree of the mesh — typically 6–8 " *
               "neighbours in 2D, so A holds ~7n non-zeros instead of n².")
    add_bullet("Reordering nodes (Cuthill–McKee, METIS) does not change " *
               "the physics; it relabels the indices and can dramatically " *
               "reduce fill-in for sparse direct solvers.")
    add_bullet("Boundary conditions are applied by editing the rows and " *
               "columns of A whose indices match the boundary node IDs — " *
               "the same node-to-row mapping is reused.")
    add_bullet("Higher-order elements (P2, P3, …) add edge/face degrees of " *
               "freedom; the connectivity table just gets more rows, but " *
               "the global-DOF-to-matrix-row rule is unchanged.")

    add_page_break()

    # ─── STEP 6: PACKAGE DESIGN PATTERNS ────────────────────────────────
    add_heading("Step 8 — Package Design Patterns: Structs, Dispatch & Plots (was Step 6)", 1)
    add_para("A well-designed Julia package separates three concerns into three " *
             "types: the grid (where), the problem (what), and the solution " *
             "(result). Multiple dispatch then lets you write one function name " *
             "— solve, plot — and Julia picks the right method automatically " *
             "based on which types you pass in. This section builds a minimal " *
             "PDESolvers package from scratch to show that workflow.")

    add_heading("6.1  The Three-Type Architecture", 2)
    add_para("Grid resolution, physical parameters, and computed results change " *
             "at different rates and for different reasons. Separating them into " *
             "their own types makes it trivial to refine the mesh, swap the " *
             "boundary condition, or compare two physics models without rewriting " *
             "the solver logic:")
    add_table(["Type",        "Holds",                         "Changes when"],
              [["Grid1D",     "n, dx, x array",                "you want finer resolution"],
               ["HeatProblem","α, T, u0, a Grid1D",            "you change the physics or IC"],
               ["Solution1D", "x, u, t, the problem reference","solve() returns one of these"]])

    add_heading("Beginner Concept Check: Structs, Instances, and Functions", 3)
    add_para("If you are new to custom data types, the interaction between structs, " *
             "variables, and functions can feel abstract. Let's clarify exactly " *
             "how they work together:")
    add_para("1. What is a Struct? A struct is a custom composite data type that " *
             "bundles multiple pieces of related data (called fields) together under " *
             "a single name. Instead of passing ten individual variables into every " *
             "function, we pack them into one structured container blueprint.")
    add_para("2. How to create a Struct? We write the struct keyword followed by " *
             "the container name, list its internal fields with their types, and end " *
             "with end. Using the Base.@kwdef macro above it automatically writes " *
             "a keyword constructor with convenient default values:")
    add_code(["Base.@kwdef struct SimpleHeatConfig",
              "    L::Float64 = 1.0       # Domain length",
              "    α::Float64 = 0.01      # Thermal diffusivity",
              "    N::Int     = 100       # Grid resolution",
              "end"])
    add_para("3. What is the variable that is the Struct (an Instance)? The code " *
             "above is just a blueprint. To use it, we call the struct's name like " *
             "a function to create a concrete object in computer memory, called an " *
             "instance, and assign it to a variable. We can then read its internal " *
             "fields using dot notation:")
    add_code(["# Create an instance; omitted fields automatically use defaults",
              "my_config = SimpleHeatConfig(α=0.05, N=200)",
              "",
              "# Access individual bundled fields using a dot (.):",
              "println(\"Configured diffusivity: \", my_config.α)",
              "println(\"Configured grid points: \", my_config.N)"])
    add_para("4. How to create a Function to receive the Struct? To write a function " *
             "that operates on our configuration, we define a function that accepts " *
             "a single argument annotated with our custom struct type. Inside the " *
             "function body, we use dot notation to unpack the exact fields we need:")
    add_code(["function print_stability_limit(config::SimpleHeatConfig)",
              "    dx = config.L / (config.N + 1)",
              "    max_dt = 0.5 * dx^2 / config.α",
              "    println(\"Maximum stable time step dt: \", max_dt)",
              "end",
              "",
              "# Pass our struct variable into the function",
              "print_stability_limit(my_config)"])
    add_para("5. Structs vs. Functions: Nouns vs. Verbs. A struct holds passive state " *
             "and data (e.g., domain dimensions, arrays, physical constants). A " *
             "function defines active behavior and logic (e.g., loops, time-stepping, " *
             "mathematical updates). Structs are the materials; functions are the " *
             "machines that transform them.")
    add_para("6. Global Scope vs. Local Scope (Crucial for Performance). Variables " *
             "defined directly in the main script body live in the global scope. Because " *
             "global variables can change their data type at any time, the Julia " *
             "compiler cannot optimize them, making loops in global scope very slow. " *
             "Variables defined inside a function body live in a local scope. The compiler " *
             "knows their exact types, compiling local operations to extremely fast " *
             "machine code. Best Practice: Always wrap your solver loops and computations " *
             "inside functions!")
    add_para("7. Where do you run these scripts? Once you write your code in a file " *
             "(e.g., solver.jl), you can run it in three ways: (a) Terminal execution: " *
             "run 'julia solver.jl' directly from your Linux command line. (b) REPL " *
             "inclusion: open the Julia REPL and type 'include(\"solver.jl\")' to load " *
             "and execute the file. (c) VS Code execution: open the file in VS Code and " *
             "press Shift+Enter to run specific blocks interactively.")

    add_heading("6.2  The Grid Struct", 2)
    add_para("Grid1D stores everything about the spatial domain. Computed fields " *
             "(dx, x) are derived from n and L so callers cannot accidentally " *
             "create an inconsistent grid. Base.@kwdef generates a keyword " *
             "constructor with defaults — Grid1D() gives a 200-point unit " *
             "interval and Grid1D(n=800, L=2.0) overrides only what you need:")
    add_code(["# src/grids.jl",
              "Base.@kwdef struct Grid1D",
              "    n  :: Int             = 200",
              "    L  :: Float64         = 1.0",
              "    dx :: Float64         = L / (n + 1)",
              "    x  :: Vector{Float64} = collect(range(dx, L - dx, length=n))",
              "end"])

    add_heading("6.3  Problem Structs", 2)
    add_para("Each physics model gets its own struct. Both embed a Grid1D so " *
             "the solver can always reach the mesh through the problem — no " *
             "extra arguments needed. Greek-letter field names (α, c) match " *
             "the mathematical notation used in the PDE:")
    add_code(["# src/problems.jl",
              "Base.@kwdef struct HeatProblem",
              "    grid :: Grid1D   = Grid1D()",
              "    α    :: Float64  = 0.01              # thermal diffusivity",
              "    T    :: Float64  = 1.0               # final time",
              "    u0   :: Function = x -> exp(-100*(x - 0.5)^2)",
              "end",
              "",
              "Base.@kwdef struct WaveProblem",
              "    grid :: Grid1D   = Grid1D()",
              "    c    :: Float64  = 1.0               # wave speed",
              "    T    :: Float64  = 1.5",
              "    u0   :: Function = x -> sin.(π .* x)",
              "end"])
    add_page_break()

    add_heading("6.4  The Solution Type and solve() Dispatch", 2)
    add_para("Solution1D is parameterised by the problem type P so that plot() " *
             "can dispatch on it at compile time — zero runtime cost. Two " *
             "solve() methods share the same name; Julia selects the right one " *
             "from the argument type alone:")
    add_code(["# src/solvers.jl",
              "struct Solution1D{P}",
              "    prob :: P",
              "    x    :: Vector{Float64}",
              "    u    :: Vector{Float64}",
              "    t    :: Float64",
              "end",
              "",
              "function solve(prob::HeatProblem)",
              "    g  = prob.grid",
              "    dt = 0.4 * g.dx^2 / prob.α          # CFL: r = α·dt/dx² ≤ 0.5",
              "    u  = prob.u0.(g.x)",
              "    for _ in 1:ceil(Int, prob.T / dt)",
              "        u[2:end-1] .+= (prob.α * dt / g.dx^2) .*",
              "            (u[3:end] .- 2 .* u[2:end-1] .+ u[1:end-2])",
              "    end",
              "    return Solution1D(prob, g.x, u, prob.T)",
              "end",
              "",
              "function solve(prob::WaveProblem)",
              "    g  = prob.grid",
              "    dt = g.dx / prob.c                  # CFL = 1 (exact for leap-frog)",
              "    λ  = (prob.c * dt / g.dx)^2",
              "    u_prev = prob.u0.(g.x)",
              "    u_curr = copy(u_prev); u_next = similar(u_prev)",
              "    for _ in 1:ceil(Int, prob.T / dt)",
              "        u_next[2:end-1] .= 2 .* u_curr[2:end-1] .- u_prev[2:end-1] .+",
              "            λ .* (u_curr[3:end] .- 2 .* u_curr[2:end-1] .+ u_curr[1:end-2])",
              "        u_next[[1,end]] .= 0.0          # Dirichlet BCs",
              "        u_prev, u_curr, u_next = u_curr, u_next, u_prev",
              "    end",
              "    return Solution1D(prob, g.x, u_curr, prob.T)",
              "end"])

    add_heading("6.5  Plot Dispatch in the Project Style", 2)
    add_para("Extending Plots.plot for your own types uses the same dispatch. " *
             "Import the function first so you add a method to the existing " *
             "generic rather than shadow it. The kwargs... passthrough lets " *
             "callers override any keyword — color, linestyle, legend — without " *
             "touching the defaults:")
    add_code(["# src/visualization.jl",
              "import Plots: plot",
              "",
              "function plot(sol::Solution1D{HeatProblem}; kwargs...)",
              "    plot(sol.x, sol.u;",
              "        xlabel = \"x\",",
              "        ylabel = \"u(x,t)\",",
              "        title  = \"Heat equation  α=$(sol.prob.α)  t=$(sol.t)\",",
              "        lw     = 2,",
              "        label  = \"Numerical\",",
              "        kwargs...)",
              "end",
              "",
              "function plot(sol::Solution1D{WaveProblem}; kwargs...)",
              "    plot(sol.x, sol.u;",
              "        xlabel = \"x\",",
              "        ylabel = \"u(x,t)\",",
              "        title  = \"Wave equation  c=$(sol.prob.c)  t=$(sol.t)\",",
              "        lw     = 2,",
              "        label  = \"Numerical\",",
              "        kwargs...)",
              "end"])
    add_page_break()

    add_heading("6.6  The Full Package Layout", 2)
    add_para("Splitting the four concerns into four files keeps each file under " *
             "200 lines and lets you include only the solver files you need. " *
             "The module entry point is the only file that contains using and " *
             "export — sub-files just define functions and types:")
    add_code(["PDESolvers/",
              "├── Project.toml",
              "└── src/",
              "    ├── PDESolvers.jl      # module entry — using, include, export",
              "    ├── grids.jl           # Grid1D",
              "    ├── problems.jl        # HeatProblem, WaveProblem",
              "    ├── solvers.jl         # solve(::HeatProblem), solve(::WaveProblem)",
              "    └── visualization.jl   # plot(::Solution1D{P})"])
    add_code(["# src/PDESolvers.jl",
              "module PDESolvers",
              "",
              "using Plots",
              "",
              "include(\"grids.jl\")",
              "include(\"problems.jl\")",
              "include(\"solvers.jl\")",
              "include(\"visualization.jl\")",
              "",
              "export Grid1D, HeatProblem, WaveProblem, Solution1D, solve",
              "",
              "end # module PDESolvers"])

    add_heading("6.7  Calling the Package", 2)
    add_para("With the types and dispatch in place, switching between physics " *
             "models is a one-word change at the call site — the grid, the " *
             "plot call, and the output pipeline stay identical. This is the " *
             "payoff of the three-type architecture:")
    add_code(["using PDESolvers, Plots",
              "",
              "# -- Heat: 400-point grid, low diffusivity, run to T=0.5 --",
              "grid  = Grid1D(n=400, L=2.0)",
              "hprob = HeatProblem(grid=grid, α=0.005, T=0.5)",
              "hsol  = solve(hprob)          # dispatches to HeatProblem method",
              "plot(hsol)                    # dispatches to Solution1D{HeatProblem}",
              "savefig(\"heat_solution.png\")",
              "",
              "# -- Wave: same grid, different physics, identical API --",
              "wprob = WaveProblem(grid=grid, c=1.5, T=2.0)",
              "wsol  = solve(wprob)          # dispatches to WaveProblem method",
              "plot(wsol)                    # dispatches to Solution1D{WaveProblem}",
              "savefig(\"wave_solution.png\")"])

    add_page_break()

    # ─── STEP 8: HPC SSH ────────────────────────────────────────────────
    add_heading("Step 9 — Connect to the HPC Cluster", 1)
    add_code(["ssh-keygen -t ed25519 -C \"your_email@example.com\"",
              "ssh-copy-id user@cluster.example.edu",
              "ssh user@cluster.example.edu",
              "module load julia/1.10",
              "julia --version"])
    add_para("Sync your project (Project.toml + Manifest.toml) to the cluster:")
    add_code(["rsync -avz ~/julia-pde/ user@cluster:~/julia-pde/",
              "ssh user@cluster",
              "cd ~/julia-pde && julia --project=. -e 'using Pkg; Pkg.instantiate()'"])
    add_page_break()

    # ─── STEP 9: SLURM ──────────────────────────────────────────────────
    add_heading("Step 10 — Submitting Julia PDE Jobs with Slurm", 1)
    add_table(["Command", "Action"],
              [["sbatch job.sh",    "Submit job"],
               ["squeue -u \$USER", "Show your jobs"],
               ["scancel JOBID",    "Cancel a job"],
               ["sinfo",            "Show partitions"],
               ["sacct -j JOBID",   "Detailed stats after run"]])
    add_para("Example heat equation job script:")
    add_code(["#!/bin/bash",
              "#SBATCH --job-name=heat_eq",
              "#SBATCH --cpus-per-task=4",
              "#SBATCH --mem=4G",
              "#SBATCH --time=00:30:00",
              "#SBATCH --output=heat_%j.out",
              "",
              "module load julia/1.10",
              "julia --threads=4 --project=. heat_equation.jl"])

    add_page_break()
    add_heading("Troubleshooting Quick Reference", 1)
    add_table(["Problem", "Solution"],
              [["wsl --install fails",     "Run as Admin; restart Windows"],
               ["julia: command not found","source ~/.bashrc"],
               ["Pkg precompile is slow",  "Normal first time; subsequent loads are fast"],
               ["SSH publickey denied",    "Re-run ssh-copy-id; ask instructor"],
               ["Slurm OOM kill",          "Increase --mem or reduce N"],
               ["CFL instability (NaN)",   "Reduce dt; check stability condition"]])

    add_page_break()
    # ─── STEP 11: GIT & GITHUB ──────────────────────────────────────────
    add_heading("Step 11 — Share Your Code with Git & GitHub", 1)
    add_para("You now have a working Julia project with PDE solvers and a " *
             "Slurm job script. The last step is to put your code under " *
             "version control with Git, and push a copy to GitHub so you can " *
             "access it from anywhere (including the HPC cluster) and share " *
             "it with classmates or your advisor.")
    add_para("What is Git? Git records snapshots of your project over time " *
             "so you can see what changed, undo mistakes, and collaborate " *
             "without emailing zip files. GitHub is a website that hosts " *
             "those snapshots online.")

    add_heading("11.1  Install Git (inside WSL)", 2)
    add_para("Open your Ubuntu terminal and run:")
    add_code(["sudo apt update",
              "sudo apt install -y git",
              "git --version   # should print: git version 2.43.x or similar"])

    add_heading("11.2  Tell Git who you are (one-time setup)", 2)
    add_para("Every snapshot is stamped with your name and email. Set them " *
             "once and Git remembers them for every project on this computer.")
    add_code(["git config --global user.name  \"Your Name\"",
              "git config --global user.email \"you@example.com\""])
    add_para("Use the same email as your GitHub account when you create one " *
             "in step 11.6.")

    add_heading("11.3  Initialise the repo in your project folder", 2)
    add_para("Go to the folder that holds your Julia project (the one with " *
             "Project.toml) and turn it into a Git repository:")
    add_code(["cd ~/projects/julia-docx-tutorial   # adjust to your folder",
              "git init",
              "git status"])
    add_para("git init creates a hidden .git/ folder where Git stores all " *
             "history. git status shows every file Git can see; right now " *
             "they will all be listed as untracked.")

    add_heading("11.4  Tell Git which files to ignore", 2)
    add_para("You do NOT want to commit generated PNG figures, the large " *
             ".docx, Julia's package cache, or editor temp files. Create a " *
             "file named .gitignore at the top of the project with this " *
             "content (one pattern per line):")
    add_code(["# Julia",
              "Manifest.toml",
              "*.jl.cov",
              "*.jl.mem",
              ".julia/",
              "",
              "# Generated artifacts",
              "output/",
              "site/output/figures/",
              "*.docx",
              "*.png",
              "",
              "# OS / editor",
              ".DS_Store",
              ".vscode/",
              "*.swp"])
    add_para("Anything matching one of these patterns is hidden from Git. " *
             "You can always remove a line later if you decide you do want " *
             "to track that file.")

    add_heading("11.5  Make your first commit", 2)
    add_para("A commit is a labelled snapshot of every file currently staged. " *
             "Stage all files, then commit with a short message:")
    add_code(["git add .",
              "git status                          # confirm what is staged",
              "git commit -m \"Initial project: Julia PDE solvers and setup guide\""])
    add_para("Git replies with something like '[main (root-commit) 9e4c42c] " *
             "Initial project ...'. That hash is the snapshot's unique ID; " *
             "you can always come back to it later.")
    add_para("Good commit messages describe what you did and why, in the " *
             "present tense. 'Add heat-equation solver and convergence test' " *
             "is better than 'changes' or 'wip'.")

    add_heading("11.6  Create a GitHub account and an empty repository", 2)
    add_bullet("Go to https://github.com and sign up (free).")
    add_bullet("Click + → New repository in the top-right corner.")
    add_bullet("Name it (for example julia-pde-tutorial), keep it public or private as you prefer.")
    add_bullet("Do NOT tick 'Add a README' or 'Add .gitignore' — your local folder already has those.")
    add_bullet("Click Create repository. GitHub shows a page with the URL of your new (empty) repo.")

    add_heading("11.7  Connect your local repo to GitHub and push", 2)
    add_para("On the GitHub repo page, copy the HTTPS URL (looks like " *
             "https://github.com/<you>/julia-pde-tutorial.git). Then in your " *
             "terminal:")
    add_code(["git branch -M main",
              "git remote add origin https://github.com/<your-username>/julia-pde-tutorial.git",
              "git push -u origin main"])
    add_para("The first time you push, GitHub will ask for your username and " *
             "a personal access token (PAT) instead of a password. Create one " *
             "at https://github.com/settings/tokens → Generate new token " *
             "(classic), give it the 'repo' scope, and paste it where Git " *
             "asks for the password. Git will cache it so you don't have to " *
             "enter it again.")

    add_heading("11.8  Day-to-day loop", 2)
    add_para("From now on, every time you change something you want to keep:")
    add_code(["git status                          # what changed?",
              "git add path/to/file.jl             # stage specific files (or 'git add .' for all)",
              "git commit -m \"Fix CFL bug in wave solver\"",
              "git push                            # send the new commits to GitHub"])
    add_para("Refresh your GitHub repo page in the browser — you should see " *
             "every file from your project, your commit message at the top, " *
             "and a clickable history. From here you can 'git clone' the " *
             "same repo on the HPC cluster and run the exact same code there.")

    # ─── SAVE ────────────────────────────────────────────────────────────
    out_path = joinpath(CFG.output_dir, CFG.docx_name)
    doc.save(out_path)
    @info "DOCX saved: $out_path"
    return out_path
end

# =============================================================================
# MAIN  — generate figures, then docx
# =============================================================================

function main()
    @info "Output directory: $(CFG.output_dir)"
    generate_all_figures()
    out = build_docx()
    @info "✓ All done. Open: $out"
end

# Only run the full pipeline when invoked as a script (`julia generate_lecture_doc.jl`);
# `include`-ing this file (e.g. to reuse a single fig_* function) won't trigger it.
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
