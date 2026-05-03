Base.@kwdef struct WaveEquation
    nx::Int = 300
    nt::Int = 300
    c::Float64 = 1.0
    L::Float64 = 1.0
    T::Float64 = 1.5
    bc::Symbol = :Dirichlet  # :Dirichlet, :Neumann, :Periodic, or :Absorbing
    f_init::Function = x -> sin.(pi * x) # Initial condition function
    dx::Float64 = L / (nx - 1)
    dt::Float64 = T / nt
    x::Any = range(0, L, length=nx)
end

WaveEquation(nx, nt, c, L, T) = WaveEquation(nx=nx, nt=nt, c=c, L=L, T=T)

function apply_boundary_conditions!(u_next, u_curr, problem::WaveEquation; is_initial=false)
    if problem.bc == :Dirichlet
        u_next[1] = 0.0
        u_next[end] = 0.0
    elseif problem.bc == :Neumann
        u_next[1] = u_next[2]
        u_next[end] = u_next[end-1]
    elseif problem.bc == :Periodic
        u_next[1] = u_next[end-1]
        u_next[end] = u_next[2]
    elseif problem.bc == :Absorbing && !is_initial
        cfl = (problem.c * problem.dt) / problem.dx
        u_next[1] = u_curr[2] + ((cfl - 1) / (cfl + 1)) * (u_next[2] - u_curr[1])
        u_next[end] = u_curr[end-1] + ((cfl - 1) / (cfl + 1)) * (u_next[end-1] - u_curr[end])
    end
end

function solve_wave_1d(problem::WaveEquation; return_history::Bool=false)
    dx = problem.dx
    dt = problem.dt
    x = problem.x

    # Step 5 — CFL:  λ = c·dt/dx ≤ 1
    λ = (problem.c * dt / dx)^2       # λ² used in update rule

    # Step 1 — initial condition: evaluate f_init over the grid
    u_curr = problem.f_init.(x)

    apply_boundary_conditions!(u_curr, u_curr, problem, is_initial=true)

    # Step 5 — zero initial velocity → u_prev = u_curr
    u_prev = copy(u_curr)
    u_next = copy(u_curr) # Copy to initialize boundary elements properly

    if return_history
        U = zeros(length(x), problem.nt + 1)
        U[:, 1] = u_curr
    end

    for i in 1:problem.nt
        # Steps 2+3+4 — leap-frog:
        #   u^{n+1} = 2u^n - u^{n-1} + λ²·(stencil)
        u_next[2:end-1] .= 2 .* u_curr[2:end-1] .- u_prev[2:end-1] .+ λ .* (u_curr[3:end] .- 2 .* u_curr[2:end-1] .+ u_curr[1:end-2])

        # Apply Boundary Conditions
        apply_boundary_conditions!(u_next, u_curr, problem)

        u_prev, u_curr = copy(u_curr), copy(u_next)

        if return_history
            U[:, i+1] = u_curr
        end
    end

    if return_history
        return collect(x), U
    else
        return collect(x), u_curr
    end
end

function plot_wave_1d(problem::WaveEquation; exact::Bool=false)
    x, u = solve_wave_1d(problem)
    p = plot(x, u, label="Numerical", xlabel="x", ylabel="u(x,T)", title="Wave equation  T=$(problem.T)", lw=2)

    if exact
        x_exact, u_exact = wave_1d_exact(problem)
        plot!(p, x_exact, u_exact, label="Exact", ls=:dash, lw=2)
    end
    return p
end

function animate_wave_1d(problem::WaveEquation; exact::Bool=false, fps::Int=30, filename="wave_animation.gif", skip::Int=1)
    x, U = solve_wave_1d(problem, return_history=true)

    # Calculate global y-limits so the plot window doesn't bounce around
    ymin, ymax = minimum(U), maximum(U)
    dy = max(0.1, (ymax - ymin) * 0.1)

    # Use skip parameter to avoid generating thousands of plots for large nt
    anim = @animate for i in 1:skip:size(U, 2)
        t = (i - 1) * problem.dt
        p = plot(x, U[:, i], ylims=(ymin - dy, ymax + dy), xlabel="x", ylabel="u",
            title="Wave Equation t=$(round(t, digits=2))", label="Numerical", lw=2)

        if exact
            prob_t = WaveEquation(nx=problem.nx, nt=problem.nt, c=problem.c, L=problem.L, T=t)
            # Override bc and f_init directly since the positional constructor might not pass them
            prob_t = WaveEquation(nx=problem.nx, nt=problem.nt, c=problem.c, L=problem.L, T=t, bc=problem.bc, f_init=problem.f_init)
            x_exact, u_exact = wave_1d_exact(prob_t)
            plot!(p, x_exact, u_exact, label="Exact", ls=:dash, lw=2)
        end
    end

    return gif(anim, filename, fps=fps)
end

function wave_1d_exact(problem::WaveEquation)
    x = collect(problem.x)
    c = problem.c
    L = problem.L
    T = problem.T

    f_init = problem.f_init

    function odd_extension(xi)
        xi_mod = mod(xi, 2L)
        if xi_mod <= L
            return f_init(xi_mod)
        else
            return -f_init(2L - xi_mod)
        end
    end

    function even_extension(xi)
        xi_mod = mod(xi, 2L)
        if xi_mod <= L
            return f_init(xi_mod)
        else
            return f_init(2L - xi_mod)
        end
    end

    function periodic_extension(xi)
        xi_mod = mod(xi, L)
        return f_init(xi_mod)
    end

    function infinite_domain(xi)
        return (xi >= 0.0 && xi <= L) ? f_init(xi) : 0.0
    end

    if problem.bc == :Neumann
        u_exact = [0.5 * (even_extension(xi - c * T) + even_extension(xi + c * T)) for xi in x]
    elseif problem.bc == :Periodic
        u_exact = [0.5 * (periodic_extension(xi - c * T) + periodic_extension(xi + c * T)) for xi in x]
    elseif problem.bc == :Absorbing
        u_exact = [0.5 * (infinite_domain(xi - c * T) + infinite_domain(xi + c * T)) for xi in x]
    else
        u_exact = [0.5 * (odd_extension(xi - c * T) + odd_extension(xi + c * T)) for xi in x]
    end

    return x, u_exact
end

function plot_wave_1d_exact(problem::WaveEquation)
    x, u = wave_1d_exact(problem)
    plot(x, u, label="Exact", xlabel="x", ylabel="u(x,T)", title="Wave equation (Exact) T=$(problem.T)", lw=2)
end