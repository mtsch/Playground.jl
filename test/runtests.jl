using Playground
using Test
using LinearAlgebra

using Distances
using IterTools
using LightGraphs
using SimpleWeightedGraphs
using NearestNeighbors
using StaticArrays
using Plots
using RecipesBase

n = 1000
# Datasets:
# Two perpendicular lines at 0.5 distance apart.
const two_lines = (vcat([@SVector[0.0, rand()] for _ in 1:n÷2],
                        [@SVector[rand(), 0.5] for _ in 1:n÷2]),
                   Euclidean(),
                   0.4)
# 2-Torus
const torus = let θ = rand(n), φ = rand(n), R = 3, r = 1
    ([@SVector[(R+r*cospi(2θ[i]))cospi(2φ[i]),
               (R+r*cospi(2θ[i]))sinpi(2φ[i]),
               r*sinpi(2θ[i])] for i in 1:n],
     Euclidean(),
     1.0)
end
# 3-sphere
const sphere3 = ([normalize(@SVector randn(4)) for _ in 1:n],
                 Euclidean(),
                 1.0)
# Integer points on a grid, manhattan distance.
const grid2 = ([@SVector[i, j, k] for i in 1:10 for j in 1:10 for k in 1:10],
               Cityblock(),
               5)

const datasets = (two_lines, torus, sphere3, grid2)

include("triangleiterator.jl")
include("geodesiccomplex.jl")
