# ── 1D ────────────────────────────────────────────────────────────────────────
Base.@kwdef struct WaveEquation <: HyperbolicProblem
    nx::Int      = 300
    nt::Int      = 300
    c::Float64   = 1.0
    L::Float64   = 1.0
    T::Float64   = 1.5
    bc::Symbol   = :Dirichlet   # :Dirichlet | :Neumann | :Periodic | :Absorbing
    f_init::Function = x -> sin.(pi * x)
    dx::Float64  = L / (nx - 1)
    dt::Float64  = T / nt
end

WaveEquation(nx, nt, c, L, T) = WaveEquation(nx=nx, nt=nt, c=c, L=L, T=T)

# ── Internal helper ────────────────────────────────────────────────────────────
function _apply_bc!(u_next, u_curr, p::WaveEquation; is_initial=false)
    if p.bc == :Dirichlet
        u_next[1] = 0.0;  u_next[end] = 0.0
    elseif p.bc == :Neumann
        u_next[1] = u_next[2];  u_next[end] = u_next[end-1]
    elseif p.bc == :Periodic
        u_next[1] = u_next[end-1];  u_next[end] = u_next[2]
    elseif p.bc == :Absorbing && !is_initial
        cfl = p.c * p.dt / p.dx
        u_next[1]   = u_curr[2]   + ((cfl - 1) / (cfl + 1)) * (u_next[2]   - u_curr[1])
        u_next[end] = u_curr[end-1] + ((cfl - 1) / (cfl + 1)) * (u_next[end-1] - u_curr[end])
    end
end

# ── Solver ─────────────────────────────────────────────────────────────────────
function solve(p::WaveEquation; return_history::Bool=false)
    x = collect(range(0, p.L, length=p.nx))
    λ = (p.c * p.dt / p.dx)^2          # CFL²; leap-frog stable when λ ≤ 1

    u_curr = p.f_init.(x)
    _apply_bc!(u_curr, u_curr, p, is_initial=true)
    u_prev = copy(u_curr)
    u_next = copy(u_curr)

    if return_history
        U = zeros(p.nx, p.nt + 1)
        U[:, 1] = u_curr
    end

    for i in 1:p.nt
        u_next[2:end-1] .= 2 .* u_curr[2:end-1] .- u_prev[2:end-1] .+
            λ .* (u_curr[3:end] .- 2 .* u_curr[2:end-1] .+ u_curr[1:end-2])
        _apply_bc!(u_next, u_curr, p)
        u_prev, u_curr = copy(u_curr), copy(u_next)
        return_history && (U[:, i+1] = u_curr)
    end

    if return_history
        t_grid = collect(range(0, p.T, length=p.nt + 1))
        return PDESolution((x, t_grid), U, p.T, p)   # 2D: space × time
    else
        return PDESolution((x,), u_curr, p.T, p)      # 1D: final snapshot
    end
end

# ── Exact solution (d'Alembert) ────────────────────────────────────────────────
function wave_exact(p::WaveEquation)
    x = collect(range(0, p.L, length=p.nx))
    f = p.f_init

    odd_ext(xi)      = (xi_m = mod(xi, 2p.L); xi_m <= p.L ? f(xi_m) : -f(2p.L - xi_m))
    even_ext(xi)     = (xi_m = mod(xi, 2p.L); xi_m <= p.L ? f(xi_m) :  f(2p.L - xi_m))
    periodic_ext(xi) = f(mod(xi, p.L))
    infinite(xi)     = (0 <= xi <= p.L) ? f(xi) : 0.0

    u_exact = if p.bc == :Neumann
        [0.5 * (even_ext(xi - p.c*p.T) + even_ext(xi + p.c*p.T)) for xi in x]
    elseif p.bc == :Periodic
        [0.5 * (periodic_ext(xi - p.c*p.T) + periodic_ext(xi + p.c*p.T)) for xi in x]
    elseif p.bc == :Absorbing
        [0.5 * (infinite(xi - p.c*p.T) + infinite(xi + p.c*p.T)) for xi in x]
    else  # Dirichlet
        [0.5 * (odd_ext(xi - p.c*p.T) + odd_ext(xi + p.c*p.T)) for xi in x]
    end

    return PDESolution((x,), u_exact, p.T, p)
end

# ── Animation ──────────────────────────────────────────────────────────────────
function animate_wave(p::WaveEquation; exact::Bool=false, fps::Int=30,
                      filename="wave_animation.gif", skip::Int=1)
    sol = solve(p, return_history=true)   # sol is a 2D PDESolution (space × time)
    x, t_grid = sol.x, sol.y
    U = sol.u

    ymin, ymax = minimum(U), maximum(U)
    dy = max(0.1, (ymax - ymin) * 0.1)

    anim = @animate for i in 1:skip:size(U, 2)
        t_now = t_grid[i]
        plot(x, U[:, i], ylims=(ymin - dy, ymax + dy),
             xlabel="x", ylabel="u", lw=2, label="Numerical",
             title="Wave t=$(round(t_now, digits=2))")
        if exact
            p_t = WaveEquation(nx=p.nx, nt=p.nt, c=p.c, L=p.L, T=t_now, bc=p.bc, f_init=p.f_init)
            ex  = wave_exact(p_t)
            plot!(ex.x, ex.u, label="Exact", ls=:dash, lw=2)
        end
    end
    return gif(anim, filename, fps=fps)
end

# ── 2D ────────────────────────────────────────────────────────────────────────
Base.@kwdef struct WaveEquation2D <: HyperbolicProblem
    Nx::Int = 100
    Ny::Int = 100
    Nt::Int = 300
    c::Float64 = 1.0
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    T::Float64 = 1.0
    f_init::Function = (x, y) -> exp(-80 * ((x - 0.5)^2 + (y - 0.5)^2))
end

function solve(p::WaveEquation2D)
    dx = p.Lx / (p.Nx - 1)
    dy = p.Ly / (p.Ny - 1)
    dt = p.T / p.Nt
    λx = (p.c * dt / dx)^2
    λy = (p.c * dt / dy)^2

    x = collect(range(0, p.Lx, length=p.Nx))
    y = collect(range(0, p.Ly, length=p.Ny))

    u_prev = [p.f_init(xi, yi) for xi in x, yi in y]
    u_curr = copy(u_prev)
    u_next = similar(u_curr)

    for _ in 1:p.Nt
        u_next[2:end-1, 2:end-1] .=
            2 .* u_curr[2:end-1, 2:end-1] .- u_prev[2:end-1, 2:end-1] .+
            λx .* (u_curr[3:end,   2:end-1] .- 2 .* u_curr[2:end-1, 2:end-1] .+ u_curr[1:end-2, 2:end-1]) .+
            λy .* (u_curr[2:end-1, 3:end  ] .- 2 .* u_curr[2:end-1, 2:end-1] .+ u_curr[2:end-1, 1:end-2])
        u_next[1, :] .= 0;  u_next[end, :] .= 0
        u_next[:, 1] .= 0;  u_next[:, end] .= 0
        u_prev = u_curr;  u_curr = u_next
        u_next = similar(u_curr)
    end
    return PDESolution((x, y), u_curr, p.T, p)
end

# ── 3D ────────────────────────────────────────────────────────────────────────
Base.@kwdef struct WaveEquation3D <: HyperbolicProblem
    Nx::Int = 40
    Ny::Int = 40
    Nz::Int = 40
    Nt::Int = 100
    c::Float64 = 1.0
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 1.0
    T::Float64 = 0.5
    f_init::Function = (x, y, z) -> exp(-80 * ((x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2))
end

function solve(p::WaveEquation3D)
    dx = p.Lx / (p.Nx - 1)
    dy = p.Ly / (p.Ny - 1)
    dz = p.Lz / (p.Nz - 1)
    dt = p.T / p.Nt
    λx = (p.c * dt / dx)^2
    λy = (p.c * dt / dy)^2
    λz = (p.c * dt / dz)^2

    x = collect(range(0, p.Lx, length=p.Nx))
    y = collect(range(0, p.Ly, length=p.Ny))
    z = collect(range(0, p.Lz, length=p.Nz))

    u_prev = [p.f_init(xi, yi, zi) for xi in x, yi in y, zi in z]
    u_curr = copy(u_prev)
    u_next = similar(u_curr)

    for _ in 1:p.Nt
        u_next[2:end-1, 2:end-1, 2:end-1] .=
            2 .* u_curr[2:end-1,2:end-1,2:end-1] .- u_prev[2:end-1,2:end-1,2:end-1] .+
            λx .* (u_curr[3:end,   2:end-1, 2:end-1] .- 2 .* u_curr[2:end-1,2:end-1,2:end-1] .+ u_curr[1:end-2, 2:end-1, 2:end-1]) .+
            λy .* (u_curr[2:end-1, 3:end,   2:end-1] .- 2 .* u_curr[2:end-1,2:end-1,2:end-1] .+ u_curr[2:end-1, 1:end-2, 2:end-1]) .+
            λz .* (u_curr[2:end-1, 2:end-1, 3:end  ] .- 2 .* u_curr[2:end-1,2:end-1,2:end-1] .+ u_curr[2:end-1, 2:end-1, 1:end-2])
        u_next[1,:,:] .= 0;  u_next[end,:,:] .= 0
        u_next[:,1,:] .= 0;  u_next[:,end,:] .= 0
        u_next[:,:,1] .= 0;  u_next[:,:,end] .= 0
        u_prev = u_curr;  u_curr = u_next
        u_next = similar(u_curr)
    end
    return PDESolution((x, y, z), u_curr, p.T, p)
end
