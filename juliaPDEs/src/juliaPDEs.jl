module juliaPDEs

using Plots, Revise, Oceananigans
using SparseArrays, LinearAlgebra

include("types.jl")
include("heat.jl")
include("wave.jl")
include("Navier_Stokes.jl")
include("poisson.jl")
include("plots.jl")

export PDEProblem, ParabolicProblem, HyperbolicProblem, EllipticProblem, IncompressibleNSProblem
export PDESolution, Grid, TestGrid
export HeatEquation, WaveEquation, PoissonEquation, NavierStokes, LidCavityFlow
export solve, solve_implicit
export animate_navier_stokes
export l2_error, convergence_table
export interior_grid, endpoint_grid
export fig_heat_equation, fig_wave_equation, fig_poisson_equation, fig_navier_stokes
export plot_solution


end # module juliaPDEs
