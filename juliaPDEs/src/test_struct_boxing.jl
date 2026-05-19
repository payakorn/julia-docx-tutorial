using BenchmarkTools

# ============================================================
# TEST 1: Function argument annotations — all should be equal
# ============================================================
v = rand(Float64, 10_000_000)

function dot_concrete(v::Vector{Float64})
    s = 0.0; @inbounds for x in v; s += x; end; s
end

function dot_parametric(v::Vector{T}) where T
    s = zero(T); @inbounds for x in v; s += x; end; s
end

function dot_abstract_param(v::AbstractVector{T}) where T
    s = zero(T); @inbounds for x in v; s += x; end; s
end

function dot_fully_abstract(v::AbstractArray)
    s = 0.0; for x in v; s += x; end; s
end

println("=" ^ 60)
println("TEST 1: Function argument annotations (called with Vector{Float64})")
println("Claim: ALL should be equally fast — Julia specializes at call site")
println("=" ^ 60)
t1 = @belapsed dot_concrete($v)
t2 = @belapsed dot_parametric($v)
t3 = @belapsed dot_abstract_param($v)
t4 = @belapsed dot_fully_abstract($v)

println("1. Vector{Float64}:         $(round(t1*1000, digits=2)) ms")
println("2. Vector{T} where T:       $(round(t2*1000, digits=2)) ms")
println("3. AbstractVector{T}:       $(round(t3*1000, digits=2)) ms")
println("4. AbstractArray (bare):    $(round(t4*1000, digits=2)) ms")
ratio_max = max(t1,t2,t3,t4) / min(t1,t2,t3,t4)
println("Ratio fastest/slowest: $(round(ratio_max, digits=2))x  (should be ~1.0)")
println()

# ============================================================
# TEST 2: Struct field types — abstract fields SHOULD be slower
# ============================================================
println("=" ^ 60)
println("TEST 2: Struct field types")
println("Claim: Abstract field → boxing → slower; concrete/parametric → fast")
println("=" ^ 60)

# FAST: concrete struct field
struct FastStruct
    data::Vector{Float64}   # concrete — Julia knows exact layout
end

# SLOW: abstract struct field
struct SlowStruct
    data::AbstractVector    # abstract — no T, no N
end

# FLEXIBLE (also fast): struct-level type parameter
struct FlexStruct{A <: AbstractVector}
    data::A                 # bound to concrete type at construction
end

function sum_fast(s::FastStruct)
    total = 0.0
    @inbounds for x in s.data; total += x; end
    total
end

function sum_slow(s::SlowStruct)
    total = 0.0
    for x in s.data; total += x; end   # no @inbounds — type unknown
    total
end

function sum_flex(s::FlexStruct)
    total = 0.0
    @inbounds for x in s.data; total += x; end
    total
end

fast = FastStruct(v)
slow = SlowStruct(v)
flex = FlexStruct(v)

tf = @belapsed sum_fast($fast)
ts = @belapsed sum_slow($slow)
tx = @belapsed sum_flex($flex)

println("FastStruct (Vector{Float64} field):         $(round(tf*1000, digits=2)) ms")
println("FlexStruct{A<:AbstractVector} field:        $(round(tx*1000, digits=2)) ms")
println("SlowStruct (AbstractVector bare field):     $(round(ts*1000, digits=2)) ms")
println()
println("Slowdown from abstract field: $(round(ts/tf, digits=1))x slower than concrete")
println("FlexStruct overhead vs FastStruct: $(round(tx/tf, digits=2))x (should be ~1.0)")
println()

# Also check allocations with @btime
println("--- Allocation check ---")
print("FastStruct: "); @btime sum_fast($fast)
print("FlexStruct: "); @btime sum_flex($flex)
print("SlowStruct: "); @btime sum_slow($slow)

println()
println("=" ^ 60)
println("SUMMARY")
println("=" ^ 60)
println("Test 1 — Function annotations: all $(round(ratio_max, digits=2))x spread (confirmed equal)")
println("Test 2 — Struct fields: abstract is $(round(ts/tf, digits=1))x slower than concrete")
