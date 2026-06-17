module juliaPDEs

using Plots, Revise, Oceananigans
using SparseArrays, LinearAlgebra
using CairoMakie
import GLMakie          # qualified-only import: provides the GLMakie window backend
                        # without clashing with CairoMakie's re-exported names
using DelimitedFiles

include("types.jl")
include("history.jl")
include("heat.jl")
include("wave.jl")
include("Navier_Stokes.jl")
include("poisson.jl")
include("plots.jl")
include("makie.jl")

export PDEProblem, ParabolicProblem, HyperbolicProblem, EllipticProblem, IncompressibleNSProblem
export PDESolution, Grid, TestGrid
export HeatEquation, WaveEquation, PoissonEquation, NavierStokes, LidCavityFlow
export solve, solve_implicit
export animate_navier_stokes
export l2_error, convergence_table
export interior_grid, endpoint_grid
export fig_heat_equation, fig_wave_equation, fig_poisson_equation, fig_navier_stokes
export plot_solution
export testmakie, makie_simple_line, makie_sin_cos_tan, makie_scatter_line
export makie_subplots, makie_contour, makie_contourf, makie_advanced
export save_sharp, save_makie_examples, sharp_display, interactive
export save_solution, load_history, SolutionWriter, save_step!, write_meta!, default_run_name


end # module juliaPDEs
