module juliaPDEs

using Plots, Revise

include("heat.jl")
include("wave.jl")

export solve_heat_1d, solve_wave_1d

end # module juliaPDEs
