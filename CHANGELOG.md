# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

---

## [0.3.0] — 2026-05-25

### Added
- **`src/plots.jl`** — new file adding figure functions for every PDE type:
  - `fig_heat_equation(p::HeatEquation{1})` — 2-panel: snapshots + space-time heatmap
  - `fig_wave_equation(p::WaveEquation{1})` — 2-panel: snapshots + characteristic-line heatmap
  - `fig_poisson_equation(p::PoissonEquation{2})` — 3-panel: f(x,y) / u(x,y) / 3-D surface
  - `fig_navier_stokes(p::LidCavityFlow)` — 2-panel: cavity schematic + speed heatmap with streamlines
  - `plot_solution(sol::PDESolution)` — generic dispatch on `sol.problem` type
- **`LidCavityFlow` struct** in `src/Navier_Stokes.jl` — lightweight streamfunction-vorticity
  (ψ-ω) lid-driven cavity solver; no Oceananigans dependency required
- **`Plots.plot` extensions** — `plot(prob)` and `plot(sol)` now dispatch to the correct
  multi-panel figure automatically for all four PDE types
- All new symbols exported: `LidCavityFlow`, `fig_heat_equation`, `fig_wave_equation`,
  `fig_poisson_equation`, `fig_navier_stokes`, `plot_solution`

### Changed
- Web tutorial pages updated to use the new one-liner API (`plot(prob; savepath="…")`)
  replacing the verbose inline `fig_*()` function definitions in `pde-wave.html`,
  `pde-poisson.html`, and `pde-ns.html`
- NS page description updated to note that the ψ-ω solver now lives in the package
- Site footer now shows the package version alongside the commit hash
- `CLAUDE.md` updated to document v0.3.0 architecture, all exports, and `plot()` usage

---

## [0.2.0] — 2026-05-19

### Added
- **`TestGrid`-driven `HeatEquation`** — dimension-free (1D/2D/3D), supports time-dependent
  Dirichlet BCs via `bc1`/`bc2` tuples
- **`solve_implicit(p::HeatEquation; θ)`** — Crank–Nicolson (θ=0.5) and backward-Euler
  (θ=1.0) implicit time stepping; unconditionally stable
- **Multi-page docs site** — split from single `index.html` into per-topic pages:
  `pde-heat.html`, `pde-wave.html`, `pde-poisson.html`, `pde-ns.html`, `pde-bc-verify.html`,
  `pde-workflow.html`, and supporting pages
- **Gmsh chapter** as a standalone page (`gmsh.html`) with mesh-to-matrix figure
- **`convergence_table`** helper in `src/poisson.jl` for empirical order-of-accuracy studies
- **`Grid` and `TestGrid` structs** in `src/types.jl` for declarative problem setup

### Changed
- `HeatEquation` refactored from ad-hoc parameters to `TestGrid`-based construction
- `PoissonEquation` and `WaveEquation` unified under the same `{N, F}` parametric pattern
- Tutorial section numbering fixed (7.x / 8.x dense numbering resolved)

---

## [0.1.0] — 2026-05 (initial releases)

### Added
- Core PDE solver package `juliaPDEs` with:
  - `HeatEquation` / `WaveEquation` — forward-Euler and leap-frog 1-D solvers
  - `PoissonEquation` — N-D sparse direct solve via Kronecker-sum Laplacian
  - `NavierStokes` — Oceananigans-backed lid-driven cavity
  - `PDESolution <: AbstractArray` — unified solution wrapper with `sol.x`, `sol.y`, `sol.u`
  - `l2_error` for Poisson verification
- `generate_lecture_doc.jl` — generates all PDE figures and assembles `.docx` via `python-docx`
- `serv.jl` — LiveServer-based static file server on port 9000
- `pde_explorer.jl` — Pluto.jl interactive notebook with live sliders
- Docker support with optional Cloudflare tunnel sidecar
- Single-page web tutorial (`site/index.html`) with KaTeX and Prism.js
- 2-D and 3-D heat/wave/Poisson solvers with parametric `{N, F}` structs
- Abstract type hierarchy: `PDEProblem → ParabolicProblem / HyperbolicProblem / EllipticProblem / IncompressibleNSProblem`
