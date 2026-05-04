function solve_navier_stokes()
    # Lid-driven cavity at Re = 100
    grid = Oceananigans.RectilinearGrid(size=(64, 64), x=(0, 1), y=(0, 1),
        topology=(Bounded, Bounded, Flat))

    # Lid-driven boundary condition: top lid moves at speed 1.0
    u_bcs = Oceananigans.FieldBoundaryConditions(north = Oceananigans.ValueBoundaryCondition(1.0))

    model = Oceananigans.NonhydrostaticModel(
        grid;
        advection=Oceananigans.UpwindBiased(order=5),
        timestepper=:RungeKutta3,
        boundary_conditions=(u=u_bcs,),
        closure=Oceananigans.ScalarDiffusivity(ν=1e-2))   # ν = 1/Re = 0.01

    Oceananigans.set!(model, u=0, v=0)

    simulation = Oceananigans.Simulation(model, Δt=1e-3, stop_time=10.0)

    # We will compute the speed field directly using Oceananigans
    speed_field = Oceananigans.Field(sqrt(model.velocities.u^2 + model.velocities.v^2))

    # Save history to memory for animation
    speed_hist = []
    
    function save_state(sim)
        Oceananigans.compute!(speed_field)
        push!(speed_hist, copy(Oceananigans.interior(speed_field, :, :, 1)))
    end

    # Save 10 frames per second of simulation time (every 0.1s = 100 steps)
    simulation.callbacks[:save] = Oceananigans.Callback(save_state, Oceananigans.IterationInterval(100))

    Oceananigans.run!(simulation)
    
    return (speed_hist=speed_hist, Δt=0.1) # 100 * 1e-3 = 0.1s per frame
end

function animate_navier_stokes(history; fps::Int=30, filename="navier_stokes_animation.gif", skip::Int=1)
    speed_hist = history.speed_hist
    dt = history.Δt
    
    anim = @animate for i in 1:skip:length(speed_hist)
        heatmap(speed_hist[i]',
            xlabel="x", ylabel="y",
            title="Navier-Stokes Speed t=$(round((i - 1) * dt, digits=2))",
            color=:viridis, aspect_ratio=1)
    end
    
    return gif(anim, filename, fps=fps)
end