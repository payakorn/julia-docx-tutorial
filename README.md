# Julia HPC Lecture Document Generator

A self-contained Julia script that generates the **Julia for HPC — Student Setup Guide** document, including:

- All 6 PDE figures (Heat, Wave, Poisson, Navier-Stokes, pipeline diagram, classification chart)
- The complete `.docx` document with embedded figures, tables, and code blocks

The whole document is generated from **a single Julia file** with all parameters in one `CONFIG` block at the top — change a number, re-run, get a new document.

---

## Quick Start

```bash
# 1. Clone or copy these files into a folder
ls
# generate_lecture_doc.jl
# Project.toml
# README.md

# 2. Activate the project and instantiate (first time only)
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 3. Run the generator
julia --project=. generate_lecture_doc.jl
```

The first run takes ~3 minutes (precompiling Plots, installing python-docx). Subsequent runs are ~30 seconds.

Output is saved to `output/`:
```
output/
├── julia_hpc_setup_guide.docx       <-- the document
└── figures/
    ├── pde_pipeline.png
    ├── pde_types.png
    ├── heat_equation.png
    ├── wave_equation.png
    ├── poisson_equation.png
    ├── navier_stokes.png
    └── convergence_study.png        <-- O(h²) verification
```

## Companion: Interactive Pluto Notebook

For live in-class demonstrations, use **`pde_explorer.jl`** — a Pluto.jl notebook with sliders for every parameter:

```bash
julia -e 'using Pkg; Pkg.add(["Pluto", "PlutoUI"]); using Pluto; Pluto.run()'
```

Then open `pde_explorer.jl` from the Pluto interface. The notebook covers:
- Heat equation with sliders for α, N, t, pulse width
- Wave equation with sliders for c, N, t, and choice of initial shape
- Poisson equation with sliders for grid size and source mode (m, n)
- Live convergence study showing O(h²) on log-log plot
- Discussion prompts for the class

---

## Customizing the Document

All parameters live in the `Config` struct at the top of `generate_lecture_doc.jl`. Edit the values, save, re-run.

### Heat Equation
```julia
heat_N         = 200          # grid points
heat_alpha     = 0.01         # thermal diffusivity α
heat_L         = 1.0          # domain length
heat_T_final   = 1.0          # final time
heat_cfl_safety= 0.4          # dt = safety * dx² / α (must be < 0.5)
heat_snapshots = [0.0, 0.05, 0.2, 0.5, 1.0]
heat_init      = x -> exp(-100*(x-0.5)^2)   # initial condition
```

### Wave Equation
```julia
wave_N         = 300
wave_c         = 1.0          # wave speed
wave_T_final   = 1.5
wave_cfl_safety= 0.4          # dt = safety * dx / c (must be < 1)
wave_snapshots = [0.1, 0.3, 0.5, 0.8, 1.2]
wave_init      = x -> x < 0.5 ? 2x : 2(1-x)  # triangular pluck
```

### Poisson Equation
```julia
poisson_N     = 50            # interior points per side
poisson_L     = 1.0
poisson_f     = (x,y) -> 2*π^2 * sin(π*x) * sin(π*y)   # source term
poisson_exact = (x,y) -> sin(π*x) * sin(π*y)            # known solution
```

### Navier-Stokes
```julia
ns_Re      = 100.0    # Reynolds number
ns_U_lid   = 1.0      # lid velocity
ns_grid_n  = 40       # streamline visualization grid
```

### Document metadata
```julia
doc_title    = "Julia for HPC"
doc_subtitle = "Student Setup Guide"
doc_author   = "Course Instructor"
```

---

## Example Customizations

**Use a sharper Gaussian for the heat equation:**
```julia
heat_init = x -> exp(-500*(x-0.5)^2)
```

**Switch to a sine-wave initial condition for the wave equation:**
```julia
wave_init = x -> sin(π*x)
```

**Demonstrate higher-order spatial frequencies in Poisson:**
```julia
poisson_f     = (x,y) -> 8*π^2 * sin(2π*x) * sin(2π*y)
poisson_exact = (x,y) -> sin(2π*x) * sin(2π*y)
```

**Higher Reynolds number cavity flow:**
```julia
ns_Re = 1000.0
```

---

## File Structure

```
generate_lecture_doc.jl
│
├── CONFIG                    Single struct holding all knobs
│
├── PDE SOLVERS              (also used in lecture as reference)
│    ├── solve_heat_1d
│    ├── solve_wave_1d
│    └── solve_poisson_2d
│
├── FIGURE GENERATION
│    ├── fig_pde_pipeline
│    ├── fig_pde_types
│    ├── fig_heat_equation
│    ├── fig_wave_equation
│    ├── fig_poisson_equation
│    └── fig_navier_stokes
│
└── DOCX GENERATION          via python-docx (auto-installed)
     └── build_docx
```

---

## How It Works

1. **Figures**: Pure Julia using `Plots.jl` with the GR backend in headless mode (`GKSwstype=100`) so it works on HPC nodes without a display.
2. **DOCX**: The `python-docx` library is called via `PythonCall.jl`. The script auto-installs python-docx with `pip` on first run.
3. **Embedded images**: Figures are saved as PNGs, then `doc.add_picture()` embeds each one with a caption.

---

## Requirements

- **Julia 1.9+** (1.10 LTS recommended)
- **Internet access** on first run to install packages
- **~500 MB disk** for the package depot

The script is self-bootstrapping — it installs missing packages (including python-docx) automatically.

---

## Troubleshooting

**`PythonCall: ImportError: No module named 'docx'`**
The script tries to install python-docx automatically. If that fails, run:
```bash
python3 -m pip install python-docx
```

**Figures look different from the docx-embedded version**
Plots.jl chooses fonts available on your system. To match the published document, install LaTeX or set:
```julia
default(fontfamily="DejaVu Sans")
```

**`MethodError: no method matching ...` after editing Config**
Make sure your `mutable struct Config` field types match what you assign. For example, `Float64` fields cannot accept `Int` literals without conversion.

---

## License

MIT — use freely for teaching.
