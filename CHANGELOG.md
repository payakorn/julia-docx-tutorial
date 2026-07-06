# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

---

## [0.5.0] — 2026-07-06

### Added
- **`CoupledHeatEquation`** (`src/coupled.jl`) — the first multi-field problem type: two
  linearly-coupled diffusing fields, `∂u/∂t = α₁∇²u − κ(u−v)` and `∂v/∂t = α₂∇²v − κ(v−u)`,
  solved as *one* sparse block linear system per time step (backward Euler, reusing the
  same `nd_laplacian` + `lu` machinery as `PoissonEquation` / `HeatEquation.solve_implicit`)
  instead of two independent per-field solves. Verified against
  `HeatEquation.solve_implicit` at κ=0 with a matching θ (exact match).
- **`CoupledPDESolution`** (`src/types.jl`) — lightweight two-field solution type returned
  by `solve(::CoupledHeatEquation)`.
- **`fig_coupled_heat(p::CoupledHeatEquation{1})`** (`src/plots.jl`) — 2-panel figure: both
  fields' final profiles overlaid, and their pointwise coupling residual `|u−v|`. Wired
  into the `Plots.plot` dispatch.
- New section **8.9 Generalizing to Coupled Systems** (`site/polish.html`, mirrored in the
  `.docx` guide as Step 8, §6.8) — why independent per-field time-stepping can't correctly
  capture a coupling term shared by both equations, and how the same one-block-solve
  technique connects to `LidCavityFlow`'s ψ/ω splitting.

### Changed
- Bumped package version to `0.5.0` in `Project.toml`.

---

## [0.4.0] — 2026-06-17

### Added
- **Solution run history** — opt-in snapshot saving wired into `solve` for `HeatEquation`
  and `WaveEquation` via `save_every` / `save_dir` (the fast path is unchanged when they
  are omitted).
- **Language-agnostic JSON format** — each run lands in its own descriptive, timestamped
  folder. *All* metadata (problem, params, grid, dt, shape, frame index) lives once in
  `meta.json`; the per-iteration files hold only the flat, column-major field array, so
  they stay tiny and load anywhere (e.g. `np.array(json.load(f)).reshape(shape, order="F")`).
- **`src/history.jl`** helpers — `save_solution`, `load_history`, and the lower-level
  `SolutionWriter` / `save_step!` / `write_meta!` / `default_run_name`.
- New docs page `solutions-io.html` — plotting a `PDESolution` and saving / reloading run
  history.

### Changed
- Added `JSON` and `Dates` dependencies to the package.

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
