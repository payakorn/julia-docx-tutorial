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
const REQUIRED_PKGS = ["Plots", "PythonCall"]

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

# Set headless backend for HPC / CI
ENV["GKSwstype"] = "100"
gr()

# =============================================================================
# CONFIG  — Edit this block to change any parameter
# =============================================================================
Base.@kwdef mutable struct Config
    # ---- output paths ----
    output_dir   :: String = "output"
    figures_dir  :: String = "output/figures"
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

function fig_navier_stokes()
    @info "Generating Figure 6: Navier-Stokes (lid-driven cavity)"

    # Left panel: cavity schematic
    p1 = plot(legend=false, framestyle=:none, aspect_ratio=:equal,
              xlim=(-0.2, 1.2), ylim=(-0.2, 1.4),
              title="Lid-Driven Cavity (NS benchmark)",
              titlefontcolor=parse(Colorant, BLUE))
    # walls
    plot!(p1, [0,1,1,0,0], [0,0,1,1,0], color=:black, lw=2.5)
    # lid arrow
    plot!(p1, [0.15, 0.85], [1.05, 1.05], arrow=true, color=:firebrick, lw=2.5)
    annotate!(p1, 0.5, 1.18,
              text("moving lid  U = $(CFG.ns_U_lid)", :firebrick, :center, 11, :bold))
    # wall labels
    annotate!(p1, -0.07, 0.5, text("u=v=0", parse(Colorant, BLUE), :right, 10))
    annotate!(p1, 1.07, 0.5,  text("u=v=0", parse(Colorant, BLUE), :left, 10))
    annotate!(p1, 0.5, -0.10, text("u=v=0", parse(Colorant, BLUE), :center, 10))
    # vortex circles
    θ = range(0, 2π, length=100)
    for r in [0.15, 0.25, 0.35]
        plot!(p1, 0.5 .+ r .* cos.(θ), 0.55 .+ r .* sin.(θ),
              color=:seagreen, lw=1, alpha=0.6)
    end
    annotate!(p1, 0.5, 0.55,
              text("vortex", :seagreen, :center, 11, :bold))

    # Right panel: streamlines (illustrative analytical approximation)
    n = CFG.ns_grid_n
    xg = range(0, 1, length=n)
    yg = range(0, 1, length=n)
    X = [x for y in yg, x in xg]
    Y = [y for y in yg, x in xg]
    # crude stream function for visualization
    Ψ = @. X*(1-X)*Y^2*(1-Y) * sin(π*X)
    # contour plot of stream function = streamlines
    p2 = contour(xg, yg, Ψ, levels=20, c=:plasma, fill=false, lw=1.2,
                 xlabel="x", ylabel="y",
                 title="Streamlines (illustrative)",
                 titlefontcolor=parse(Colorant, BLUE),
                 aspect_ratio=:equal, colorbar_title="Ψ")

    plt = plot(p1, p2, layout=(1,2),
               size=(CFG.fig_width, CFG.fig_height+50), dpi=CFG.fig_dpi,
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

function generate_all_figures()
    @info "=== Generating all figures ==="
    fig_pde_pipeline()
    fig_pde_types()
    fig_heat_equation()
    fig_wave_equation()
    fig_poisson_equation()
    fig_navier_stokes()
    fig_convergence_study()
    @info "All figures saved to $(CFG.figures_dir)"
end

# =============================================================================
# DOCX GENERATION  — uses python-docx via PythonCall
# =============================================================================

"""
Loads python-docx (auto-installs via pip if missing) and returns the module.
"""
function load_python_docx()
    # Try to import; if it fails, install via pip
    try
        return pyimport("docx")
    catch
        @info "python-docx not found — installing via pip..."
        python_exe = pyconvert(String, pyimport("sys").executable)
        run(`$python_exe -m pip install --quiet python-docx`)
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

    # ─── OVERVIEW ───────────────────────────────────────────────────────
    add_heading("Overview & Prerequisites", 1)
    add_para("This guide walks you through setting up a complete Julia " *
             "environment for scientific computing and HPC on Windows. " *
             "Steps 1–3 install Julia. Step 4 introduces PDE packages. " *
             "Step 5 walks the full PDE-solving workflow with worked " *
             "examples for the heat, wave, Poisson, and Navier-Stokes " *
             "equations. Steps 6–8 connect to the HPC cluster.")
    add_heading("Time Estimate", 2)
    add_table(["Step", "Estimated Time"],
              [["1 — Install WSL2",          "10–15 min"],
               ["2 — juliaup + Julia",       "5–10 min"],
               ["3 — Install Julia Packages","10–15 min"],
               ["4 — Learn PDE Packages",    "15–20 min"],
               ["5 — PDE Workflow + Examples","30–60 min"],
               ["6 — Configure VS Code",     "5–10 min"],
               ["7 — HPC SSH Setup",         "10 min"],
               ["8 — Submit Slurm Jobs",     "10–15 min"],
               ["Total",                     "~95–155 min"]])
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

    # ─── STEP 3: PACKAGES ───────────────────────────────────────────────
    add_heading("Step 3 — Installing Julia Packages", 1)
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
    add_heading("Step 4 — Key Julia Packages for PDE Solving", 1)
    add_table(["Package", "Use case"],
              [["LinearAlgebra (stdlib)",     "Dense matrix ops"],
               ["SparseArrays (stdlib)",      "Sparse PDE matrices"],
               ["FFTW.jl",                    "Spectral methods"],
               ["DifferentialEquations.jl",   "ODE/PDE solver suite"],
               ["MethodOfLines.jl",           "Auto-discretization of PDEs"],
               ["KrylovKit.jl",               "GMRES, CG iterative solvers"],
               ["Ferrite.jl / Gridap.jl",     "Finite element methods"],
               ["Oceananigans.jl",            "Production CFD on HPC"],
               ["Plots.jl / Makie.jl",        "Visualization"],
               ["MPI.jl / CUDA.jl",           "HPC parallelism"]])
    add_page_break()

    # ─── STEP 5: PDE WORKFLOW ───────────────────────────────────────────
    add_heading("Step 5 — The PDE Solving Workflow", 1)
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

    # ─── STEP 6: VS CODE ────────────────────────────────────────────────
    add_heading("Step 6 — Configure VS Code for Julia + WSL", 1)
    add_para("Install these VS Code extensions:")
    add_table(["Extension", "Purpose"],
              [["WSL",            "Connect VS Code to WSL"],
               ["Julia",          "Syntax, IntelliSense, REPL"],
               ["Jupyter",        "Run notebooks with Julia kernel"],
               ["Remote — SSH",   "Connect to HPC cluster"]])
    add_page_break()

    # ─── STEP 7: HPC SSH ────────────────────────────────────────────────
    add_heading("Step 7 — Connect to the HPC Cluster", 1)
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

    # ─── STEP 8: SLURM ──────────────────────────────────────────────────
    add_heading("Step 8 — Submitting Julia PDE Jobs with Slurm", 1)
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

main()
