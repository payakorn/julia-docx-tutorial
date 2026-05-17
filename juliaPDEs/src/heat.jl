# ── 1D ────────────────────────────────────────────────────────────────────────
Base.@kwdef struct HeatEquation <: ParabolicProblem
    N::Int       = 200
    Nt::Int      = 1000            # number of time steps (dt = T / Nt)
    α::Float64   = 0.01
    L::Float64   = 1.0
    T::Float64   = 1.0
    f_init::Function = x -> exp(-100 * (x - 0.5)^2)
end

function solve(p::HeatEquation)
    dx = p.L / (p.N + 1)
    x  = collect(range(dx, p.L - dx, length=p.N))
    dt = p.T / p.Nt
    r  = p.α * dt / dx^2
    r ≤ 0.5 || @warn "HeatEquation stability: r = α·dt/dx² = $r > 0.5 — increase Nt or reduce N."

    u = p.f_init.(x)
    for _ in 1:p.Nt
        u[2:end-1] .+= r .* (u[3:end] .- 2 .* u[2:end-1] .+ u[1:end-2])
    end
    return PDESolution((x,), u, p.T, p)
end

# ── 2D ────────────────────────────────────────────────────────────────────────
Base.@kwdef struct HeatEquation2D <: ParabolicProblem
    Nx::Int      = 50
    Ny::Int      = 50
    Nt::Int      = 100             # number of time steps (dt = T / Nt)
    α::Float64   = 0.01
    Lx::Float64  = 1.0
    Ly::Float64  = 1.0
    T::Float64   = 0.5
    f_init::Function = (x, y) -> exp(-50 * ((x - 0.5)^2 + (y - 0.5)^2))
end

function solve(p::HeatEquation2D)
    dx = p.Lx / (p.Nx + 1)
    dy = p.Ly / (p.Ny + 1)
    x  = collect(range(dx, p.Lx - dx, length=p.Nx))
    y  = collect(range(dy, p.Ly - dy, length=p.Ny))
    dt = p.T / p.Nt
    rx = p.α * dt / dx^2
    ry = p.α * dt / dy^2
    rx + ry ≤ 0.5 || @warn "HeatEquation2D stability: rx + ry = $(rx + ry) > 0.5 — increase Nt."

    u     = [p.f_init(xi, yi) for xi in x, yi in y]
    u_new = similar(u)
    for _ in 1:p.Nt
        u_new[2:end-1, 2:end-1] .=
            u[2:end-1, 2:end-1] .+
            rx .* (u[3:end,   2:end-1] .- 2 .* u[2:end-1, 2:end-1] .+ u[1:end-2, 2:end-1]) .+
            ry .* (u[2:end-1, 3:end  ] .- 2 .* u[2:end-1, 2:end-1] .+ u[2:end-1, 1:end-2])
        u_new[1, :] .= 0;  u_new[end, :] .= 0
        u_new[:, 1] .= 0;  u_new[:, end] .= 0
        u .= u_new
    end
    return PDESolution((x, y), u, p.T, p)
end

# ── 3D ────────────────────────────────────────────────────────────────────────
Base.@kwdef struct HeatEquation3D <: ParabolicProblem
    Nx::Int      = 20
    Ny::Int      = 20
    Nz::Int      = 20
    Nt::Int      = 50              # number of time steps (dt = T / Nt)
    α::Float64   = 0.01
    Lx::Float64  = 1.0
    Ly::Float64  = 1.0
    Lz::Float64  = 1.0
    T::Float64   = 0.2
    f_init::Function = (x, y, z) -> exp(-50 * ((x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2))
end

function solve(p::HeatEquation3D)
    dx = p.Lx / (p.Nx + 1)
    dy = p.Ly / (p.Ny + 1)
    dz = p.Lz / (p.Nz + 1)
    x  = collect(range(dx, p.Lx - dx, length=p.Nx))
    y  = collect(range(dy, p.Ly - dy, length=p.Ny))
    z  = collect(range(dz, p.Lz - dz, length=p.Nz))
    dt = p.T / p.Nt
    rx = p.α * dt / dx^2
    ry = p.α * dt / dy^2
    rz = p.α * dt / dz^2
    rx + ry + rz ≤ 0.5 || @warn "HeatEquation3D stability: rx + ry + rz = $(rx+ry+rz) > 0.5 — increase Nt."

    u     = [p.f_init(xi, yi, zi) for xi in x, yi in y, zi in z]
    u_new = similar(u)
    for _ in 1:p.Nt
        u_new[2:end-1, 2:end-1, 2:end-1] .=
            u[2:end-1, 2:end-1, 2:end-1] .+
            rx .* (u[3:end,   2:end-1, 2:end-1] .- 2 .* u[2:end-1, 2:end-1, 2:end-1] .+ u[1:end-2, 2:end-1, 2:end-1]) .+
            ry .* (u[2:end-1, 3:end,   2:end-1] .- 2 .* u[2:end-1, 2:end-1, 2:end-1] .+ u[2:end-1, 1:end-2, 2:end-1]) .+
            rz .* (u[2:end-1, 2:end-1, 3:end  ] .- 2 .* u[2:end-1, 2:end-1, 2:end-1] .+ u[2:end-1, 2:end-1, 1:end-2])
        u_new[1, :, :] .= 0;  u_new[end, :, :] .= 0
        u_new[:, 1, :] .= 0;  u_new[:, end, :] .= 0
        u_new[:, :, 1] .= 0;  u_new[:, :, end] .= 0
        u .= u_new
    end
    return PDESolution((x, y, z), u, p.T, p)
end
