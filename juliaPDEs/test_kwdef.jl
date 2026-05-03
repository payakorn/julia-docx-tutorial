Base.@kwdef struct WaveEquation
    N::Int = 300
    c::Float64 = 1.0
    L::Float64 = 1.0
    T::Float64 = 1.5
    dx::Float64 = L / (N + 1)
    x::Any = range(dx, L - dx, length=N)
    dt::Float64 = 0.4 * dx / c
end

w = WaveEquation()
println(w.dx)
println(w.dt)
println(w.x)
