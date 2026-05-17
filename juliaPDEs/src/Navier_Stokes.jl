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
    ν = 1.0 / p.Re

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
    simulation = Oceananigans.Simulation(model, Δt=p.T / p.Nt, stop_time=p.T)

    speed_field = Oceananigans.Field(
        sqrt(model.velocities.u^2 + model.velocities.v^2))

    Oceananigans.compute!(speed_field)
    Oceananigans.run!(simulation)
    Oceananigans.compute!(speed_field)

    x = collect(range(0, 1, length=p.nx))
    y = collect(range(0, 1, length=p.ny))
    u = copy(Oceananigans.interior(speed_field, :, :, 1))

    return PDESolution((x, y), u, p.T, p)
end

# ── Animation (runs its own solve with history snapshots) ──────────────────────
function animate_navier_stokes(p::NavierStokes; fps::Int=30,
                               filename="navier_stokes_animation.gif", skip::Int=1)
    ν = 1.0 / p.Re

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
    simulation = Oceananigans.Simulation(model, Δt=p.T / p.Nt, stop_time=p.T)

    speed_field = Oceananigans.Field(
        sqrt(model.velocities.u^2 + model.velocities.v^2))

    speed_hist = Matrix{Float64}[]
    function save_frame(sim)
        Oceananigans.compute!(speed_field)
        push!(speed_hist, copy(Oceananigans.interior(speed_field, :, :, 1)))
    end

    simulation.callbacks[:save] = Oceananigans.Callback(
        save_frame, Oceananigans.IterationInterval(p.save_interval))

    Oceananigans.run!(simulation)

    Δt_frame = p.save_interval * (p.T / p.Nt)
    anim = @animate for i in 1:skip:length(speed_hist)
        heatmap(speed_hist[i]',
            xlabel="x", ylabel="y", color=:viridis, aspect_ratio=1,
            title="Navier-Stokes speed  t=$(round((i-1)*Δt_frame, digits=2))")
    end
    return gif(anim, filename, fps=fps)
end
