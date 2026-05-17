# ── Wave equation — dimension-free ────────────────────────────────────────────
#
#   ∂²u/∂t² = c² Δu        on the Cartesian box  ∏ᵢ [aᵢ, bᵢ]    in ℝᴺ
#   u = 0  on ∂Ω
#   u(x, 0) = f_init(x...),   ∂u/∂t(x, 0) = 0
#
# Same N, F parametrisation as HeatEquation: N is the spatial dimension and
# F is the initial-condition closure type.
#
Base.@kwdef struct WaveEquation{N, F} <: HyperbolicProblem
    N_grid::NTuple{N, Int} = (300,)                              # grid points per axis (defaults to 1D)
    a::NTuple{N, Float64} = ntuple(_ -> 0.0, length(N_grid))     # lower bound per axis
    b::NTuple{N, Float64} = ntuple(_ -> 1.0, length(N_grid))     # upper bound per axis
    Nt::Int      = 300                                           # number of time steps (dt = T / Nt)
    T::Float64   = 1.0                                           # final simulation time
    c::Float64   = 1.0                                           # wave speed
    f_init::F    = x -> sin(π * x)                               # default is 1D — override for higher N
end

function solve(p::WaveEquation{N, F}) where {N, F}
    # 1. Build the per-axis endpoint grids via the shared helper.
    grids = ntuple(i -> endpoint_grid(p.a[i], p.b[i], p.N_grid[i]), N)
    d           = ntuple(i -> grids[i][1], N)
    axes_coords = ntuple(i -> grids[i][2], N)

    # 2. Time step from user-supplied Nt, with a CFL check (Σ λᵢ ≤ 1 for stability).
    dt = p.T / p.Nt
    λ  = ntuple(i -> (p.c * dt / d[i])^2, N)
    sum(λ) ≤ 1.0 || @warn "WaveEquation CFL: Σλᵢ = $(sum(λ)) > 1 — increase Nt."
    nsteps = p.Nt

    # 3. Initialize three buffers (∂u/∂t = 0 at t=0 ⇒ u_prev = u_curr).
    u_curr = zeros(Float64, p.N_grid)
    for I in CartesianIndices(u_curr)
        coords = ntuple(dim -> axes_coords[dim][I[dim]], N)
        u_curr[I] = p.f_init(coords...)
    end

    # 4. Enforce Dirichlet zero boundaries ONCE at initialization.
    for i in 1:N
        selectdim(u_curr, i, 1) .= 0.0
        selectdim(u_curr, i, p.N_grid[i]) .= 0.0
    end

    u_prev = copy(u_curr)
    u_next = similar(u_curr)

    inner_range = CartesianIndices(ntuple(i -> 2:(p.N_grid[i]-1), N))
    e = ntuple(i -> CartesianIndex(ntuple(j -> j == i ? 1 : 0, N)), N)

    # 5. Main time-stepping loop — leap-frog (2nd-order in time) + (2N+1)-point Laplacian.
    for _ in 1:nsteps
        for I in inner_range
            laplacian = 0.0
            for i in 1:N
                @inbounds laplacian += λ[i] * (u_curr[I+e[i]] - 2.0 * u_curr[I] + u_curr[I-e[i]])
            end
            @inbounds u_next[I] = 2.0 * u_curr[I] - u_prev[I] + laplacian
        end
        for i in 1:N
            selectdim(u_next, i, 1) .= 0.0
            selectdim(u_next, i, p.N_grid[i]) .= 0.0
        end
        # O(1) three-way pointer rotation — the old u_prev buffer is reused as the next u_next.
        u_prev, u_curr, u_next = u_curr, u_next, u_prev
    end

    return PDESolution(axes_coords, u_curr, p.T, p)
end
