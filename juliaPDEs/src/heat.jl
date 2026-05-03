function solve_heat_1d(; N=200, α=0.01, L=1.0, T=1.0)
    # Step 1 — grid
    dx     = L / (N + 1)
    x      = range(dx, L - dx, length=N)

    # Step 5 — stable Δt:  r = α·dt/dx² ≤ 0.5
    dt     = 0.4 * dx^2 / α
    nsteps = ceil(Int, T / dt)

    # Step 1 — initial condition u(x,0)
    u = @. exp(-100 * (x - 0.5)^2)

    for _ in 1:nsteps
        # Steps 3+4 — Rule 1 (stencil) inside Rule 2 (Euler step):
        #   u_i^{n+1} = u_i^n + r·(u_{i+1} - 2u_i + u_{i-1})
        u[2:end-1] .+= (α * dt / dx^2) .*
            (u[3:end] .- 2 .* u[2:end-1] .+ u[1:end-2])
        # Dirichlet BCs: u[1] = u[end] = 0 — never updated
    end
    return collect(x), u
end


"""
# How to use this function
using Plots
x, u = solve_heat_1d()
plot(x, u, xlabel="x", ylabel="u(x,T)", title="Heat equation  T=1.0", lw=2)
"""