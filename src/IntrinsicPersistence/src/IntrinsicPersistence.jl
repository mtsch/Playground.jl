module IntrinsicPersistence

using PersistenceBarcodes

using Distances
using LightGraphs, SimpleWeightedGraphs
using LightGraphs: DijkstraState
using NearestNeighbors
using StaticArrays
using RecipesBase

include("TriangleIterator.jl")
include("result_types.jl")
include("persistence.jl")

export persistence

end
