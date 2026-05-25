# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A Julia project that generates a `.docx` teaching document ("Julia for HPC — Student Setup Guide") with embedded PDE solver figures and an interactive web interface. It solves 4 PDE types (heat, wave, Poisson, Navier-Stokes), produces PNG figures, and assembles them into a Word document via `python-docx`.

## Commands

### Local Development

```bash
# Install dependencies (first time)
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Generate the .docx and all figures
julia --project=. generate_lecture_doc.jl

# Launch web server on port 9000
julia --project=. serv.jl 9000
# Open http://localhost:9000

# Interactive PDE notebook (requires Pluto)
julia -e 'using Pluto; Pluto.run()'
# Then open pde_explorer.jl in the Pluto UI
```

### Docker

```bash
# Build and start (server on port 9000, optional cloudflare tunnel)
docker-compose up --build

# With Cloudflare tunnel
TUNNEL_TOKEN=<token> docker-compose up --build
```

## Architecture

### Key Files

| File | Role |
|------|------|
| `generate_lecture_doc.jl` | Main generator (~2073 lines) — solves PDEs, renders figures, builds `.docx` |
| `serv.jl` | `LiveServer.jl` web server serving `site/` on port 9000 |
| `pde_explorer.jl` | Pluto.jl interactive notebook with live sliders for all PDEs |
| `juliaPDEs/src/` | Reusable PDE solver package (heat, wave, Poisson, Navier-Stokes) |
| `site/index.html` | Static web interface with embedded KaTeX and Prism.js |
| `CondaPkg.toml` | Declares `python-docx` (auto-installed via pip by `PythonCall.jl`) |

### Data Flow

1. `generate_lecture_doc.jl` runs the numerical solvers and writes PNGs to `site/output/figures/`
2. It calls `python-docx` via `PythonCall.jl` to write `output/julia_hpc_setup_guide.docx`
3. `serv.jl` serves `site/` statically; `site/output/` is a symlink to `output/`

### Configuration

All tunable parameters live in the `Config` struct near the top of `generate_lecture_doc.jl` (lines 44–92): grid sizes, time steps, CFL numbers, PDE coefficients (α, c, Re), document metadata, and figure styling. Edit there to change solver behavior or document appearance.

### `juliaPDEs` Package (v0.3.0)

A companion package inside `juliaPDEs/`. Source files:

| File | Role |
|------|------|
| `src/types.jl` | Abstract type hierarchy, `PDESolution`, `Grid`, `TestGrid`, grid helpers |
| `src/heat.jl` | `HeatEquation` — forward Euler + Crank–Nicolson (`solve`, `solve_implicit`) |
| `src/wave.jl` | `WaveEquation` — leap-frog, dimension-free (1D/2D/3D) |
| `src/poisson.jl` | `PoissonEquation` — sparse direct solve, `l2_error`, `convergence_table` |
| `src/Navier_Stokes.jl` | `LidCavityFlow` (ψ-ω solver) + `NavierStokes` (Oceananigans-backed) |
| `src/plots.jl` | `plot()` dispatch for all Problem/Solution types; `fig_*` figure functions |

Key exports:
- **Structs**: `HeatEquation`, `WaveEquation`, `PoissonEquation`, `LidCavityFlow`, `NavierStokes`
- **Solvers**: `solve`, `solve_implicit`, `animate_navier_stokes`
- **Figures**: `plot(prob)` / `plot(sol)` — dispatches to the correct multi-panel figure automatically
- **Errors**: `l2_error`, `convergence_table`
- **Helpers**: `interior_grid`, `endpoint_grid`

#### `plot()` dispatch

Extending `Plots.plot` means users just call `plot(prob)` or `plot(sol)`:

```julia
using juliaPDEs, Plots

plot(HeatEquation(testgrid=tg, Nt=2000, T=1.0, α=0.01))   # 2-panel: snapshots + space-time
plot(WaveEquation(N_grid=(300,), Nt=600, T=1.5, f_init=…)) # 2-panel: snapshots + space-time
plot(PoissonEquation())                                      # 3-panel: f / u / 3D surface
plot(LidCavityFlow(Re=100.0))                               # 2-panel: schematic + speed heatmap
plot(solve(prob))                                            # same dispatch via PDESolution
```

Pass `savepath="file.png"` to save to disk. The explicit `fig_*` functions are also exported for
direct use (`fig_heat_equation`, `fig_wave_equation`, `fig_poisson_equation`, `fig_navier_stokes`).

Has its own `Project.toml` with heavier deps (CUDA, Makie, Oceananigans) not needed by the main project.

## Dependencies

- **Julia 1.10 LTS** — required
- Runtime Julia packages: `Plots` (GR backend), `PythonCall`, `Gmsh`, `LinearAlgebra`, `SparseArrays`, `LiveServer`, `PlutoUI`
- Python: `python-docx` — auto-installed on first run via `PythonCall`/CondaPkg
- Docker: `julia:1.10` base image; optional `cloudflared` sidecar
- First run takes ~3 minutes (package precompilation + pip install); subsequent runs ~30 seconds
