module juliaPDEs

using Plots, Revise, Oceananigans
using SparseArrays, LinearAlgebra

include("types.jl")
include("heat.jl")
include("wave.jl")
include("Navier_Stokes.jl")
include("poisson.jl")

export PDEProblem, ParabolicProblem, HyperbolicProblem, EllipticProblem, IncompressibleNSProblem
export PDESolution
export HeatEquationND, HeatEquation, HeatEquation2D, HeatEquation3D
export WaveEquation, WaveEquation2D, WaveEquation3D
export PoissonEquation, PoissonEquation3D
export NavierStokes
export solve
export animate_wave, animate_navier_stokes
export l2_error, convergence_table, wave_exact


end # module juliaPDEs
