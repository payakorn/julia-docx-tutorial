using BenchmarkTools

# 1. Define the struct as a subtype of AbstractVector
struct EveryOtherArray{T, A<:AbstractVector{T}} <: AbstractVector{T}
    parent::A
end

# 2. Obligatory methods for the AbstractArray interface
Base.size(v::EveryOtherArray) = (div(length(v.parent) + 1, 2),)

# When someone asks for index 'i', we actually fetch '2i - 1' from the parent
Base.getindex(v::EveryOtherArray, i::Int) = v.parent[2i - 1]

Base.IndexStyle(::Type{<:EveryOtherArray}) = IndexLinear()

# Setup data (1 million random floats)
const data = rand(1_000_000)
const custom_subset = EveryOtherArray(data)

# Benchmark standard built-in sum on a standard vector
println("Standard Vector Sum:")
@btime sum($data)

# Benchmark the exact same built-in sum on our custom struct
println("\nCustom Subset Struct Sum:")
@btime sum($custom_subset)
