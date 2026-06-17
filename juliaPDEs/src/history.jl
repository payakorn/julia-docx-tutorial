# ── Solution history — language-agnostic JSON snapshots ───────────────────────
#
# A run is saved to its own folder named descriptively + timestamped, e.g.
#
#   runs/heatequation_2d_N50x50_2026-06-17T08-30-12/
#     meta.json           ← problem, params, grid, dt, shape + index of frames
#     step_00000.json     ← field values ONLY (flat, column-major) — no metadata
#     step_00010.json
#     ...
#
# All per-snapshot metadata (iteration, time, filename) lives once in meta.json
# so the individual frame files stay as small as possible. The frame files are
# plain JSON number arrays, so any language can load a frame and reshape it with
# `meta["shape"]` (column-major / Fortran order):
#
#   Python:  np.array(json.load(open(f))).reshape(meta["shape"], order="F")
#
using Dates
import JSON

# ── Run-folder naming ─────────────────────────────────────────────────────────

"""
    default_run_name(problem) -> String

Descriptive, sortable, timestamped folder name, e.g.
`heatequation_2d_N50x50_2026-06-17T08-30-12`. Falls back gracefully for
problems that don't expose an `N_grid` field.
"""
function default_run_name(problem)
    label = lowercase(string(nameof(typeof(problem))))
    stamp = Dates.format(now(), dateformat"yyyy-mm-ddTHH-MM-SS")
    if hasproperty(problem, :N_grid)
        ng = getfield(problem, :N_grid)
        return string(label, "_", length(ng), "d_N", join(ng, "x"), "_", stamp)
    else
        return string(label, "_", stamp)
    end
end

# Axis label for meta.json: x/y/z for the usual cases, axisN beyond 3-D.
_axis_name(i::Integer) = i <= 3 ? ("x", "y", "z")[i] : "axis$i"

# ── Writer — accumulates a frame index, then writes meta.json ─────────────────

"""
    SolutionWriter(dir; problem, grid, dt, shape)

Create (mkpath) a run folder `dir` and prepare to stream field snapshots into
it. Use [`save_step!`](@ref) per frame, then [`write_meta!`](@ref) once at the
end. The `grid` is an `NTuple` of coordinate vectors (one per axis).
"""
mutable struct SolutionWriter
    dir       :: String
    problem   :: Any
    grid      :: Tuple
    dt        :: Float64
    shape     :: Tuple
    snapshots :: Vector{Dict{String,Any}}
end

function SolutionWriter(dir::AbstractString; problem, grid, dt, shape)
    mkpath(dir)
    return SolutionWriter(String(dir), problem, Tuple(grid), Float64(dt),
                          Tuple(shape), Dict{String,Any}[])
end

"""
    save_step!(w::SolutionWriter, u, step, iter, t)

Write the field `u` for frame `step` as a flat (column-major) JSON number array
to `step_<step>.json`, and record its `iter`/`t`/`file` in the frame index.
"""
function save_step!(w::SolutionWriter, u::AbstractArray, step::Integer,
                    iter::Integer, t::Real)
    fname = @sprintf("step_%05d.json", step)
    open(joinpath(w.dir, fname), "w") do io
        JSON.print(io, vec(u))            # column-major flatten; reshape via meta["shape"]
    end
    push!(w.snapshots, Dict{String,Any}(
        "step" => Int(step), "iteration" => Int(iter),
        "t" => Float64(t), "file" => fname))
    return fname
end

"""
    write_meta!(w::SolutionWriter; params = Dict())

Write `meta.json` for the run: dimension, grid coordinates, field shape, dt,
the full frame index, and any solver `params` worth recording. Call once after
all `save_step!`s.
"""
function write_meta!(w::SolutionWriter; params::AbstractDict = Dict{String,Any}())
    meta = Dict{String,Any}(
        "problem"   => string(nameof(typeof(w.problem))),
        "dimension" => length(w.grid),
        "shape"     => collect(w.shape),
        "order"     => "column-major",     # how to reshape the flat frame arrays
        "dt"        => w.dt,
        "created"   => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
        "grid"      => Dict(_axis_name(i) => collect(w.grid[i]) for i in 1:length(w.grid)),
        "params"    => Dict{String,Any}(params),
        "snapshots" => w.snapshots,
    )
    open(joinpath(w.dir, "meta.json"), "w") do io
        JSON.print(io, meta, 2)            # pretty-printed, 2-space indent
    end
    return joinpath(w.dir, "meta.json")
end

# ── Convenience: save a single (final) solution ───────────────────────────────

"""
    save_solution(sol::PDESolution, dir = joinpath("runs", default_run_name(sol.problem)))

Save one solution as a one-frame run folder (`meta.json` + `step_00000.json`).
Returns the folder path.
"""
function save_solution(sol::PDESolution,
                       dir::AbstractString = joinpath("runs", default_run_name(sol.problem)))
    p = getfield(sol, :problem)
    dt = hasproperty(p, :T) && hasproperty(p, :Nt) ? p.T / p.Nt : 0.0
    w  = SolutionWriter(dir; problem=p, grid=getfield(sol, :grid),
                        dt=dt, shape=size(sol.u))
    save_step!(w, sol.u, 0, 0, sol.t)
    write_meta!(w; params=_problem_params(p))
    return dir
end

# Pull the scalar/numeric fields of a problem struct into a JSON-friendly Dict.
function _problem_params(p)
    d = Dict{String,Any}()
    for name in propertynames(p)
        v = getfield(p, name)
        v isa Union{Real,AbstractString,Tuple} && (d[string(name)] = v isa Tuple ? collect(v) : v)
    end
    return d
end

# ── Loading back into Julia ───────────────────────────────────────────────────

"""
    load_history(dir) -> (meta, times, frames)

Read a saved run folder back into Julia. `meta` is the parsed `meta.json`,
`times` is a `Vector{Float64}`, and `frames` is a `Vector{Array{Float64}}` each
reshaped to `meta["shape"]` (column-major). Works for any language's writer as
long as it follows the same layout.
"""
function load_history(dir::AbstractString)
    meta  = JSON.parsefile(joinpath(dir, "meta.json"))
    shape = Tuple(Int.(meta["shape"]))
    snaps = meta["snapshots"]
    times  = Float64[s["t"] for s in snaps]
    frames = map(snaps) do s
        flat = Float64.(JSON.parsefile(joinpath(dir, s["file"])))
        reshape(flat, shape)
    end
    return meta, times, frames
end
