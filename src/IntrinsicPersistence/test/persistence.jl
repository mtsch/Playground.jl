@testset "persistence helpers" begin
    n = 100
    nrepeats = 10

    @testset "sampledensity" begin
        using IntrinsicPersistence.sampledensity

        for dim in 2:nrepeats+1
            # Points
            pts = rand(dim, n)
            M   = pairwise(Euclidean(), pts)
            M2  = copy(M)
            for i in 1:n
                M2[i, i] = Inf
            end
            @test @inferred(sampledensity(M)) == maximum(minimum(M2, 1))

            # Graphs
            G = Graph(n, 2n)
            @test @inferred(sampledensity(G)) == 1
            A = sprand(n, n, 0.1)
            A += A'
            G = SimpleWeightedGraph(A)
            @test @inferred(sampledensity(G)) == maximum(A)
        end
    end

    @testset "edgetoorder" begin
        using IntrinsicPersistence.edgetoorder

        for _ in 1:nrepeats
            idxs = Tuple{Int, Int}[]
            A = spzeros(n, n)
            for i in 1:n
                j, k = 0, 0
                while j == k || A[j, k] ≠ 0
                    j, k = rand(1:n, 2)
                end
                A[j, k] = i
                A[k, j] = i
                push!(idxs, (j, k))
            end

            # Make sure this appears as a unit test.
            @test @inferred(edgetoorder(SimpleWeightedGraph(A))) ≠ nothing
            ord = edgetoorder(SimpleWeightedGraph(A))
            for i in 1:n
                @test ord[idxs[i]...] == i
            end
        end
    end

    @testset "spannergraph" begin
        using IntrinsicPersistence.spannergraph

        # Metrics that don't need a parameter.
        metrics = [Distances.Chebyshev(),
                   Distances.Cityblock(),
                   Distances.Euclidean(),
                   Distances.Hamming(),
                   Distances.HellingerDist(),
                   Distances.Jaccard(),
                   Distances.MeanAbsDeviation(),
                   Distances.RMSDeviation()]

        for dim in 2:nrepeats+1, m in metrics
            pts1 = rand(dim, n)
            M1   = pairwise(m, pts1)
            pts2 = [rand(SVector{dim, Float64}) for _ in 1:n]
            M2   = pairwise(m, reinterpret(Float64, pts2, (dim, n)))

            @test @inferred(spannergraph(pts1, m, 2maximum(M1))).weights ≈ M1
            @test @inferred(spannergraph(pts2, m, 2maximum(M2))).weights ≈ M2

            r1 = 2 * sampledensity(M1)
            r2 = 2 * sampledensity(M2)
            for i in 1:n, j in 1:n
                if M1[i, j] ≥ r1
                    M1[i, j] = 0
                    M1[j, i] = 0
                end
                if M2[i, j] ≥ r2
                    M2[i, j] = 0
                    M2[j, i] = 0
                end
            end
            @test @inferred(spannergraph(pts1, m)).weights ≈ M1
            @test @inferred(spannergraph(pts2, m)).weights ≈ M2
        end
    end
end

@testset "persistence" begin
    @testset "graphs" begin

        function squarelattice(n, m)
            G = Graph(n*m)
            L = reshape(1:n*m, n, m)
            for j in 1:m, i in 1:n
                i < n && add_edge!(G, L[i,j], L[i+1,j])
                j < m && add_edge!(G, L[i,j], L[i,j+1])
                i < n && j < m && add_edge!(G, L[i,j], L[i+1,j+1])
            end
            G
        end

        function glue!(G::Graph, r1, r2)
            length(r1) == length(r2) || throw("Edges should have the same length!")
            for i in eachindex(r1)
                add_edge!(G, r1[i], r2[i])
            end
            G
        end

        # Torus:
        nt = 51; mt = 21
        torus = glue!(squarelattice(nt,mt), 1:nt, nt*mt-nt+1:nt*mt)
        torus = glue!(torus, 1:nt:nt*mt, nt:nt:nt*mt)

        # Klein bottle:
        nk = 41; mk = 19
        klein = glue!(squarelattice(nk,mk), 1:nk, nk*mk-nk+1:nk*mk)
        klein = glue!(klein, 1:nk:nk*mk, reverse(nk:nk:nk*mk))

        # Real projective plane:
        np = 31; mp = 31
        pplane = glue!(squarelattice(np,mp), 1:np, reverse(np*mp-np+1:np*mp))
        pplane = glue!(pplane, 1:np:np*mp, reverse(np:np:np*mp))

        results = []
        for (idx, space) in enumerate([torus, klein, pplane])
            res1 = persistence(space, use_hash = true)
            res2 = persistence(space, use_hash = false)

            @test length(res1.generators) == length(res2.generators)
            for i in 1:length(res1.generators)
                @test res1.generators[i].path     == res2.generators[i].path
                @test res1.generators[i].triangle == res2.generators[i].triangle
                @test length(res1.generators[i])  == length(res2.generators[i])
            end

            push!(results, res1)
        end

        # Torus:
        @test length(results[1].generators)   == 2
        @test length(results[1].generators[1]) ≈ mt atol = 2
        @test length(results[1].generators[2]) ≈ nt atol = 2

        # Klein bottle:
        @test length(results[2].generators)   == 2
        @test length(results[2].generators[1]) ≈ mk atol = 2
        @test length(results[2].generators[2]) ≈ nk atol = 2

        # Real projective plane:
        @test length(results[3].generators)   == 1
        @test length(results[3].generators[1]) ≈ mp atol = 2

        # TODO: make a shorter path and make sure it's found.
    end

    @testset "points" begin

        circ = [sin.(linspace(0, 2π, 1000))'; cos.(linspace(0, 2π, 1000))']
        res  = persistence(circ)
        @test length(res.generators)   == 1
        @test length(res.generators[1]) ≈ 2π atol = 0.01

        rnd = rand(3, 1000)
        @test length(persistence(rnd).generators) == 0

        # TODO:
        # * randomly sampled spaces with large R
    end
end
