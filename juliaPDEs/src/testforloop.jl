using BenchmarkTools

v = rand(Float64, 10_000_000)   # always a Vector{Float64}

# 1. Concrete: Vector{Float64}
function dot_concrete(v::Vector{Float64})
    s = 0.0
    @inbounds for x in v
        s += x
    end
    s
end

# 2. Parametric: Vector{T}  (element type is a type param)
function dot_parametric(v::Vector{T}) where T
    s = zero(T)
    @inbounds for x in v
        s += x
    end
    s
end

# 3. Abstract with type param: AbstractVector{T}
function dot_abstract_param(v::AbstractVector{T}) where T
    s = zero(T)
    @inbounds for x in v
        s += x
    end
    s
end

# 4. Fully abstract: AbstractArray (no T, no N)
function dot_fully_abstract(v::AbstractArray)
    s = 0.0
    for x in v
        s += x
    end   # no @inbounds — type unknown
    s
end

println("1. Vector{Float64}:           ", @btime dot_concrete($v))
println("2. Vector{T}:                 ", @btime dot_parametric($v))
println("3. AbstractVector{T}:         ", @btime dot_abstract_param($v))
println("4. AbstractArray (bare):      ", @btime dot_fully_abstract($v))