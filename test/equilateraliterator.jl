import Playground.EquilateralIterators.getsortededges

@testset "EquilateralIterator" begin
    distmats = map(datasets) do (pts, metric, r)
        pts = pts[1:100]
        T = eltype(eltype(pts))
        n = length(eltype(pts))
        m = length(pts)
        pairwise(metric, reshape(reinterpret(T, pts), (n, m))), r
    end

    @testset "sortededges" begin
        for (dists, _) in distmats

            sortededges = getsortededges(dists)

            @test issorted([dists[e...] for e in sortededges])
            @test all(>(e...) for e in sortededges)
            @test length(sortededges) == binomial(size(dists, 1), 2)
        end
    end

    @testset "unique, sorted, equilateral" begin
        for (dists, r) in distmats
            Δlt(t1, t2) = isless(maximum(dists[i, j] for (i, j) in IterTools.subsets(t1, 2)),
                                 maximum(dists[i, j] for (i, j) in IterTools.subsets(t2, 2)))
            triangles = collect(equilaterals(dists, r))

            @test allunique(triangles)
            # triangles by side length
            @test issorted(getindex.(triangles, 1), lt = Δlt)
            # triangles by r
            @test issorted(getindex.(triangles, 2))

            function isequilateral(t1)
                a, b, c = (dists[i, j] for (i, j) in IterTools.subsets(t1, 2))
                abs(a - b) ≤ r && abs(b - c) ≤ r && abs(a - c) ≤ r
            end

            @test all(isequilateral, getindex.(triangles, 1))
        end
    end

    @testset "counts" begin
        h = √3/2
        #   *
        #  * *
        # * * *
        pts = [0 1 2 0.5 1.5 1;
               0 0 0 h   h   2h]
        dists  = pairwise(Euclidean(), pts)
        Δs = collect(equilaterals(dists, 0.1))

        @test length(Δs) == 5

        #   * * *
        #  * * *
        # * * *
        h   = √3/2
        pts = [0 1 2 0.5 1.5 2.5 1  2  3;
               0 0 0 h   h   h   2h 2h 2h]
        dists  = pairwise(Euclidean(), pts)
        Δs = collect(equilaterals(dists, 0.1))

        @test length(Δs) == 12

        # All points the same distance apart -- all triangles are eqilateral.
        n = 100
        dists = ones(Int, n, n)
        for i in 1:n
            dists[i, i] = 0
        end
        Δs = collect(equilaterals(Symmetric(dists), 0.0))

        @test length(Δs) == binomial(n, 3)
    end
end
