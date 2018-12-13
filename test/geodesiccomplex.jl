@testset "GeodesicComplex" begin

    @testset "Construction & data structure properties" begin
        for (pts, metric, r) in datasets
            lmks, cover = GeodesicComplexes.getcover(nothing, pts, r, KDTree(pts, metric))
            @test all(1 .≤ lmks .≤ length(pts))
            @test reduce(union, cover) == Set(1:length(pts))
            @test issorted(lmks)
            @test all(eachindex(lmks)) do i
                all(evaluate(metric, pts[lmks[i]], pts[j]) ≤ r for j in cover[i])
            end
            gc = GeodesicComplex(pts, r, metric = metric)

            # TODO: indexing
            lmrks = landmarks(gc)
            othrs = nonlandmarks(gc)
            @test allunique(lmrks)
            @test allunique(othrs)
            @test length(intersect(lmrks, othrs)) == 0
            @test length(lmrks) + length(othrs) == length(points(gc))
            @test Set(union(lmrks, othrs)) == Set(points(gc))
            nlandmarks(gc) == length(lmrks)

            @test all(weight(e) ≈ evaluate(metric, lmrks[src(e)], lmrks[dst(e)])
                      for e in edges(gc))
            @test all(r ≤ weight(e) ≤ 2r
                      for e in edges(gc))
        end
    end

    @testset "LightGraphs interface" begin
        for (pts, metric, r) in datasets
            gc = GeodesicComplex(pts, r, metric = metric)
            # Run some functions to see that the interface is defined properly
            @test betweenness_centrality(gc) isa Array
            @test dijkstra_shortest_paths(gc, 1) isa LightGraphs.DijkstraState
        end
    end

    @testset "Plotting" begin
        # Count number of series.
        d = Dict{Symbol, Any}()
        for (pts, metric, r) in datasets
            gc = GeodesicComplex(pts, r, metric = metric)
            # ¯\_(ツ)_/¯
            @test length(RecipesBase.apply_recipe(d, gc)) == 2
        end
    end
end
