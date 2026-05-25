# ── Lid-Driven Cavity — streamfunction-vorticity (ψ-ω) formulation ──────────
#
#   ∂ω/∂t + (∂ψ/∂y)(∂ω/∂x) - (∂ψ/∂x)(∂ω/∂y) = ν∇²ω   (vorticity transport)
#   ω = −∇²ψ                                              (vorticity from ψ)
#   u = ∂ψ/∂y,   v = −∂ψ/∂x                              (velocity from ψ)
#
# Eliminates pressure — only scalar solves needed. The diffusion term is
# implicit (unconditionally stable); the only constraint is the convective
# CFL: U·dt/h ≤ 1. Wall vorticity uses Thom's formula.
#
Base.@kwdef struct LidCavityFlow <: IncompressibleNSProblem
    N::Int         = 48
    Re::Float64    = 100.0
    T::Float64     = 10.0
    U_lid::Float64 = 1.0
end

# Returns (ψ, ω, u_vel, v_vel, speed, h) — all N×N arrays.
function _solve_psi_omega(p::LidCavityFlow)
    N  = p.N
    L  = 1.0
    h  = L / (N - 1)
    ν  = p.U_lid * L / p.Re
    dt = 0.30 * h / p.U_lid
    nsteps = ceil(Int, p.T / dt)

    ψ = zeros(N, N)
    ω = zeros(N, N)

    n_int = N - 2
    nn    = n_int^2
    II = Int[]; JJ = Int[]; VV = Float64[]
    for jj in 1:n_int, ii in 1:n_int
        k = (jj - 1) * n_int + ii
        push!(II, k); push!(JJ, k); push!(VV, -4.0)
        ii > 1     && (push!(II, k); push!(JJ, k - 1);     push!(VV, 1.0))
        ii < n_int && (push!(II, k); push!(JJ, k + 1);     push!(VV, 1.0))
        jj > 1     && (push!(II, k); push!(JJ, k - n_int); push!(VV, 1.0))
        jj < n_int && (push!(II, k); push!(JJ, k + n_int); push!(VV, 1.0))
    end
    Lap    = sparse(II, JJ, VV, nn, nn) ./ h^2
    Inn    = sparse(I, nn, nn)
    A_diff = lu(Inn - dt * ν * Lap)
    Lap_lu = lu(Lap)
    coeff  = ν * dt / h^2

    for _ in 1:nsteps
        # Wall vorticity via Thom's formula
        @views ω[2:N-1, N] .= -2 .* ψ[2:N-1, N-1] ./ h^2 .- 2 * p.U_lid / h
        @views ω[2:N-1, 1] .= -2 .* ψ[2:N-1, 2]   ./ h^2
        @views ω[1, 2:N-1] .= -2 .* ψ[2,   2:N-1] ./ h^2
        @views ω[N, 2:N-1] .= -2 .* ψ[N-1, 2:N-1] ./ h^2

        # Explicit convection
        ψy  = (ψ[2:N-1, 3:N]   .- ψ[2:N-1, 1:N-2]) ./ (2h)
        ψx  = (ψ[3:N,   2:N-1] .- ψ[1:N-2, 2:N-1]) ./ (2h)
        ωx  = (ω[3:N,   2:N-1] .- ω[1:N-2, 2:N-1]) ./ (2h)
        ωy  = (ω[2:N-1, 3:N]   .- ω[2:N-1, 1:N-2]) ./ (2h)
        adv = ψy .* ωx .- ψx .* ωy

        # Implicit diffusion: (I − νΔt L) ω^{n+1} = ω^n − Δt·adv
        rhs_int = ω[2:N-1, 2:N-1] .- dt .* adv
        rhs_int[1,   :] .+= coeff .* ω[1,   2:N-1]
        rhs_int[end, :] .+= coeff .* ω[N,   2:N-1]
        rhs_int[:,   1] .+= coeff .* ω[2:N-1, 1]
        rhs_int[:, end] .+= coeff .* ω[2:N-1, N]
        ω[2:N-1, 2:N-1] .= reshape(A_diff \ vec(rhs_int), n_int, n_int)

        # Streamfunction Poisson: L ψ = −ω
        ψ[2:N-1, 2:N-1] .= reshape(Lap_lu \ -vec(ω[2:N-1, 2:N-1]), n_int, n_int)
    end

    u_vel = zeros(N, N)
    v_vel = zeros(N, N)
    u_vel[:, 2:N-1] .=  (ψ[:, 3:N]   .- ψ[:, 1:N-2])   ./ (2h)
    v_vel[2:N-1, :] .= -(ψ[3:N, :]   .- ψ[1:N-2, :])   ./ (2h)
    u_vel[:, N] .= p.U_lid
    return ψ, ω, u_vel, v_vel, sqrt.(u_vel.^2 .+ v_vel.^2), h
end

function solve(p::LidCavityFlow)
    _, _, _, _, speed, _ = _solve_psi_omega(p)
    _, x = endpoint_grid(1.0, p.N)
    _, y = endpoint_grid(1.0, p.N)
    return PDESolution((x, y), speed, p.T, p)
end

# ── Oceananigans-backed Navier-Stokes ────────────────────────────────────────

Base.@kwdef struct NavierStokes <: IncompressibleNSProblem
    nx::Int            = 64
    ny::Int            = 64
    Nt::Int            = 10000   # number of time steps (Δt = T / Nt)
    Re::Float64        = 100.0
    T::Float64         = 10.0
    lid_speed::Float64 = 1.0
    save_interval::Int = 100     # iterations between saved frames
end

function solve(p::NavierStokes)
    # 1. Time step from user-supplied Nt; kinematic viscosity from Reynolds number.
    Δt = p.T / p.Nt
    ν  = 1.0 / p.Re

    # 2. Build a 2-D rectilinear lid-driven-cavity grid via Oceananigans.
    grid = Oceananigans.RectilinearGrid(
        size=(p.nx, p.ny), x=(0, 1), y=(0, 1),
        topology=(Bounded, Bounded, Flat))

    # 3. Boundary conditions: u = lid_speed on the north wall, no-slip elsewhere.
    u_bcs = Oceananigans.FieldBoundaryConditions(
        north=Oceananigans.ValueBoundaryCondition(p.lid_speed))

    # 4. Assemble the nonhydrostatic model (advection + diffusion + BCs).
    model = Oceananigans.NonhydrostaticModel(
        grid;
        advection   = Oceananigans.UpwindBiased(order=5),
        timestepper = :RungeKutta3,
        boundary_conditions = (u=u_bcs,),
        closure     = Oceananigans.ScalarDiffusivity(ν=ν))

    Oceananigans.set!(model, u=0, v=0)
    simulation = Oceananigans.Simulation(model, Δt=Δt, stop_time=p.T)

    # 5. Main run — Oceananigans advances the fields to stop_time = p.T.
    speed_field = Oceananigans.Field(
        sqrt(model.velocities.u^2 + model.velocities.v^2))

    Oceananigans.compute!(speed_field)
    Oceananigans.run!(simulation)
    Oceananigans.compute!(speed_field)

    # 6. Extract the final speed field and wrap it as a PDESolution.
    _, x = endpoint_grid(1.0, p.nx)
    _, y = endpoint_grid(1.0, p.ny)
    u    = copy(Oceananigans.interior(speed_field, :, :, 1))

    return PDESolution((x, y), u, p.T, p)
end

# ── Animation (runs its own solve with history snapshots) ──────────────────────
function animate_navier_stokes(p::NavierStokes; fps::Int=30,
                               filename="navier_stokes_animation.gif", skip::Int=1)
    # 1. Same set-up as solve() — Δt, ν, grid, BCs, model.
    Δt = p.T / p.Nt
    ν  = 1.0 / p.Re

    grid = Oceananigans.RectilinearGrid(
        size=(p.nx, p.ny), x=(0, 1), y=(0, 1),
        topology=(Bounded, Bounded, Flat))

    u_bcs = Oceananigans.FieldBoundaryConditions(
        north=Oceananigans.ValueBoundaryCondition(p.lid_speed))

    model = Oceananigans.NonhydrostaticModel(
        grid;
        advection   = Oceananigans.UpwindBiased(order=5),
        timestepper = :RungeKutta3,
        boundary_conditions = (u=u_bcs,),
        closure     = Oceananigans.ScalarDiffusivity(ν=ν))

    Oceananigans.set!(model, u=0, v=0)
    simulation = Oceananigans.Simulation(model, Δt=Δt, stop_time=p.T)

    # 2. Snapshot the speed field every `save_interval` iterations via a callback.
    speed_field = Oceananigans.Field(
        sqrt(model.velocities.u^2 + model.velocities.v^2))

    speed_hist = Matrix{Float64}[]
    function save_frame(sim)
        Oceananigans.compute!(speed_field)
        push!(speed_hist, copy(Oceananigans.interior(speed_field, :, :, 1)))
    end

    simulation.callbacks[:save] = Oceananigans.Callback(
        save_frame, Oceananigans.IterationInterval(p.save_interval))

    # 3. Main run — populates speed_hist by the callback above.
    Oceananigans.run!(simulation)

    # 4. Build the animation, one frame per saved snapshot.
    Δt_frame = p.save_interval * Δt
    anim = @animate for i in 1:skip:length(speed_hist)
        heatmap(speed_hist[i]',
            xlabel="x", ylabel="y", color=:viridis, aspect_ratio=1,
            title="Navier-Stokes speed  t=$(round((i-1)*Δt_frame, digits=2))")
    end
    return gif(anim, filename, fps=fps)
end
