Base.@kwdef struct WaveEquation
    N::Int = 300
    c::Float64 = 1.0
    L::Float64 = 1.0
    T::Float64 = 1.5
    dx::Float64 = L / (N + 1)
    dt::Float64 = 0.4 * dx / c
    x::Any = range(dx, L - dx, length=N)
end

WaveEquation(N, c, L, T) = WaveEquation(N=N, c=c, L=L, T=T)

w = WaveEquation(1000, 2, 2, 2.5)
println(w.N, " ", w.dx)
