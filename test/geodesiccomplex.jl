using Distances
using LightGraphs
using SimpleWeightedGraphs
using NearestNeighbors
using StaticArrays
using Plots
using RecipesBase

using LinearAlgebra

const n = 10_000
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
const grid2 = ([@SVector[i, j] for i in 1:100 for j in 1:100],
               Cityblock(),
               5)

const datasets = (two_lines, torus, sphere3, grid2)

@testset "GeodesicComplex" begin

    @testset "Construction & data structure properties" begin
        for (pts, metric, r) in datasets
            landmarks, cover = GeodesicComplexes.getcover(pts, r, KDTree(pts, metric))
            @test all(1 .≤ landmarks .≤ length(pts))
            @test reduce(union, cover) == Set(1:length(pts))

            for gc in (GeodesicComplex(pts, r, metric = metric),
                       GeodesicComplex(pts, r, metric = metric, witness = false))

                landmarks = landmark_points(gc)
                witnesses = nonlandmark_points(gc)
                @test allunique(landmarks)
                @test allunique(witnesses)
                @test length(intersect(landmarks, witnesses)) == 0
                @test length(landmarks) + length(witnesses) == length(points(gc))
                n_landmarks(gc) == length(landmarks)

                @test all(weight(e) ≈ evaluate(metric, landmarks[src(e)], landmarks[dst(e)])
                          for e in edges(gc))
                @test all(weight(e) ≤ 2r
                          for e in edges(gc))
            end
        end
    end

    @testset "LightGraphs interface" begin
        for (pts, metric, r) in datasets
            for gc in (GeodesicComplex(pts, r, metric = metric),
                       GeodesicComplex(pts, r, metric = metric, witness = false))
                # Run some functions to see that the interface is defined properly
                @test betweenness_centrality(gc) isa Array
                @test dijkstra_shortest_paths(gc, 1) isa LightGraphs.DijkstraState
            end
        end
    end

    @testset "Plotting" begin
        # Count number of series.
        d = Dict{Symbol, Any}()
        for (pts, metric, r) in datasets
            for gc in (GeodesicComplex(pts, r, metric = metric),
                       GeodesicComplex(pts, r, metric = metric, witness = false))
                @test length(RecipesBase.apply_recipe(d, gc)) == 2
            end
        end
    end
end
