module juliaPDEs

using Plots, Revise, Oceananigans
using SparseArrays, LinearAlgebra

include("types.jl")
include("heat.jl")
include("wave.jl")
include("Navier_Stokes.jl")
include("poisson.jl")

export PDEProblem, ParabolicProblem, HyperbolicProblem, EllipticProblem, IncompressibleNSProblem
export PDESolution, Grid, TestGrid
export HeatEquation, WaveEquation, PoissonEquation, NavierStokes
export solve, solve_implicit
export animate_navier_stokes
export l2_error, convergence_table
export interior_grid, endpoint_grid


end # module juliaPDEs
