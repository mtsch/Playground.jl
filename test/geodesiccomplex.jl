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
