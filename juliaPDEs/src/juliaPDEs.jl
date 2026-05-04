module juliaPDEs

using Plots, Revise, Oceananigans

include("heat.jl")
include("wave.jl")
include("Navier_Stokes.jl")

export solve_heat_1d, solve_wave_1d, solve_navier_stokes, animate_wave_1d, animate_navier_stokes

end # module juliaPDEs
