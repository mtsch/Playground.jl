@testset "TriangleIterator" begin
    using IntrinsicPersistence.equilaterals

    n = 100
    nrepeats = 10

    @testset "ordertoedge" begin
        using IntrinsicPersistence.ordertoedge

        for dim in 2:nrepeats+1
            pts = rand(dim, n)
            M   = pairwise(Euclidean(), pts)
            ord = ordertoedge(M)

            @test issorted(M[i, j] for (i, j) in ord)
            @test all(i > j for (i, j) in ord)
        end
    end

    @testset "type stability and iterator traits" begin
        using IntrinsicPersistence.Triangle

        for T in [Float64, Float32, Int64, Int32]
            if T <: Integer
                pts = rand(T, 3, n) .÷ T(16)
            else
                pts = rand(T, 3, n)
            end
            M = pairwise(Cityblock(), pts) # preserves Ints.

            @test_inferred equilaterals(M)
            Δs = equilaterals(M)
            @test_inferred start(Δs)
            st = start(Δs)
            @test_inferred done(Δs, st)
            @test_inferred next(Δs, st)
            @test_inferred first(Δs)

            @test Base.iteratoreltype(Δs) == Base.HasEltype()
            @test eltype(Δs) == Triangle{T}
            @test typeof(first(Δs)) == eltype(Δs)
            @test Base.iteratorsize(Δs) == Base.SizeUnknown()
        end
    end

    @testset "unique, sorted and equilateral?" begin
        using IntrinsicPersistence.sampledensity

        for dim in 2:nrepeats+1
            pts = rand(dim, n)
            M   = pairwise(Euclidean(), pts)
            Δs  = collect(equilaterals(M))
            s   = sampledensity(M)

            # Sensible length.
            @test length(Δs) < binomial(n, 3)

            # Unique.
            @test Δs == unique(Δs)

            # Sorted.
            for i in 1:length(Δs)-1
                @test Δs[i].r ≤ Δs[i+1].r
            end

            # Equilateral.
            for Δ in Δs
                i, j, k = Δ.vertices
                a = M[i, j]
                b = M[j, k]
                c = M[k, i]

                @test isapprox(a, b, atol = 2s) &&
                      isapprox(b, c, atol = 2s) &&
                      isapprox(c, a, atol = 2s)
                @test max(a, b, c) == Δ.r
            end
        end
    end

    @testset "counts" begin
        h = √3/2
        #   *
        #  * *
        # * * *
        pts = [0 1 2 0.5 1.5 1;
               0 0 0 h   h   2h]
        M  = pairwise(Euclidean(), pts)
        # General position:
        M .+= 0.001 * Symmetric(rand(size(M)))
        Δs = collect(equilaterals(M, tol = 0.1))

        @test Δs == unique(Δs)
        @test length(Δs) == 5

        #   * * *
        #  * * *
        # * * *
        h   = √3/2
        pts = [0 1 2 0.5 1.5 2.5 1  2  3;
               0 0 0 h   h   h   2h 2h 2h]
        M  = pairwise(Euclidean(), pts)
        # General position:
        M .+= 0.001 * Symmetric(rand(size(M)))
        Δs = collect(equilaterals(M, tol = 0.1))

        @test Δs == unique(Δs)
        @test length(Δs) == 12

        #  * *
        # *   *
        #  * *
        pts = [1 2 2.5 2  1  0.5;
               0 0 h   2h 2h h]
        M  = pairwise(Euclidean(), pts)
        # General position:
        M .+= 0.001 * Symmetric(rand(size(M)))
        Δs = collect(equilaterals(M, tol = 0.1))
        @test Δs == unique(Δs)
        @test length(Δs) == 2

        # All points the same distance apart with noise.
        M = ones(n, n) + 0.001 * Symmetric(rand(n, n))
        for i in 1:n
            M[i, i] = 0
        end
        Δs = collect(equilaterals(M, tol = 0.1))

        @test Δs == unique(Δs)
        @test length(Δs) == binomial(n, 3)

        # All points the same distance apart without noise --
        # each triangle is returned three times.
        M = ones(n, n)
        for i in 1:n
            M[i, i] = 0
        end
        Δs = collect(equilaterals(M, tol = 0))

        @test length(Δs) == binomial(n, 3) * 3
    end
end
