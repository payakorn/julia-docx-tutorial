module juliaPDEs

using Plots, Revise, Oceananigans
using SparseArrays, LinearAlgebra
using BenchmarkTools
using Printf

include("types.jl")
include("history.jl")
include("heat.jl")
include("wave.jl")
include("Navier_Stokes.jl")
include("poisson.jl")
include("plotting.jl")

export PDEProblem, ParabolicProblem, HyperbolicProblem, EllipticProblem, IncompressibleNSProblem
export PDESolution
export HeatEquation, WaveEquation, PoissonEquation, NavierStokes
export solve
export animate_navier_stokes
export l2_error, convergence_table
export interior_grid, endpoint_grid
export plot_solution
export save_solution, load_history, SolutionWriter, save_step!, write_meta!, default_run_name


end # module juliaPDEs
