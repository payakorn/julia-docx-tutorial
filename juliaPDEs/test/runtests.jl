using Test
using juliaPDEs

@testset "juliaPDEs" begin

  @testset "Ind indexing" begin
    d = Ind((:u, :v, :w), (4, 5))
    for var in (:u, :v, :w), I in CartesianIndices((4, 5))
      flat = d[var, I.I...]
      (var2, I2) = d[flat]
      @test var2 == var
      @test I2 == I
    end
  end

  @testset "CoupledPDESystem ≡ CoupledHeatEquation (M = 2)" begin
    grid = Grid(a=(0.0,), b=(1.0,), stepsize=(0.02,))
    f1(x) = sin(pi * x)
    f2(x) = 0.5 * sin(2pi * x)

    che = CoupledHeatEquation(grid=grid, α=(0.01, 0.02), κ=0.5, f_init=(f1, f2), Nt=200, T=0.1)
    sol_che = solve(che)

    sys = CoupledPDESystem(grid=grid, vars=(:u, :v), α=(0.01, 0.02),
                            K=((0.0, 0.5), (0.5, 0.0)), f_init=(f1, f2), Nt=200, T=0.1)
    sol_sys = solve(sys)

    @test sol_che.u ≈ sol_sys.u
    @test sol_che.v ≈ sol_sys.v
  end

  @testset "CoupledPDESystem at K=0 ≡ HeatEquation.solve_implicit (θ=1)" begin
    grid = Grid(a=(0.0,), b=(1.0,), stepsize=(0.02,))
    f_init(x) = sin(pi * x)
    tg = TestGrid(grid, (t -> 0.0,), (t -> 0.0,), f_init)

    heat = HeatEquation(testgrid=tg, Nt=200, T=0.1, α=0.01)
    sol_heat = solve_implicit(heat; θ=1.0)

    sys = CoupledPDESystem(grid=grid, vars=(:u, :v), α=(0.01, 0.01),
                            f_init=(f_init, f_init), Nt=200, T=0.1)   # K defaults to zero coupling
    sol_sys = solve(sys)

    @test sol_sys.u ≈ sol_heat.u
    @test sol_sys.v ≈ sol_heat.u
  end

  @testset "HeatEquation smoke test" begin
    grid = Grid(a=(0.0,), b=(1.0,), stepsize=(0.02,))
    tg = TestGrid(grid, (t -> 0.0,), (t -> 0.0,), x -> sin(pi * x))
    p = HeatEquation(testgrid=tg, Nt=100, T=0.1, α=0.01)
    sol = solve(p)
    @test size(sol.u) == grid.numgrid
    @test all(isfinite, sol.u)
  end

  @testset "WaveEquation smoke test" begin
    p = WaveEquation(N_grid=(80,), a=(0.0,), b=(1.0,), Nt=200, T=0.2, c=1.0, f_init=x -> sin(pi * x))
    sol = solve(p)
    @test length(sol.u) == 80
    @test all(isfinite, sol.u)
  end

  @testset "PoissonEquation smoke test" begin
    p = PoissonEquation(N_grid=(20, 20))
    sol = solve(p)
    err = l2_error(sol)
    @test err !== nothing
    @test err < 0.1
  end

end
