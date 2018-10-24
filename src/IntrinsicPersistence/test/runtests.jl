using Base.Test
using IntrinsicPersistence
using Distances
using LightGraphs, SimpleWeightedGraphs
using StaticArrays

macro test_inferred(ex)
    :(@test @inferred($(esc(ex))) ≠ nothing)
end

include("TriangleIterator.jl")
include("persistence.jl")
