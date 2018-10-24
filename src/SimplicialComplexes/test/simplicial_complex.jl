@testset "SimplicialComplexes" begin
    @testset "maximalsimplices" begin
        sxcx = SimplicialComplex()
        push!(sxcx, (1, 2, 3, 4))
        push!(sxcx, (3, 4, 5))
        maxsxs = collect(values(SimplicialComplexes.maximalsimplices(sxcx)))
        @test [3, 4, 5] in maxsxs
        @test [1, 2, 3, 4] in maxsxs
        @test !([1, 2, 3] in maxsxs)
    end

    @testset "maximalcofaces, maximalfaces" begin
        maximalcofaces = SimplicialComplexes.maximalcofaces
        maximalfaces = SimplicialComplexes.maximalfaces

        sxcx = SimplicialComplex()
        push!(sxcx, (1, 2, 3, 4)) # lab = 1
        push!(sxcx, (2, 3, 4, 5)) # lab = 2
        push!(sxcx, (3, 4, 6))    # lab = 3

        @test maximalcofaces(sxcx, (1, 2)) == Set(1)
        @test maximalcofaces(sxcx, (4, 5)) == Set(2)
        @test maximalcofaces(sxcx, (3, 4)) == Set([1, 2, 3])

        @test maximalfaces(sxcx, (3, 4, 5, 6)) == Set(3)
        @test maximalfaces(sxcx, 1:10) == Set([1, 2, 3])
        @test maximalfaces(sxcx, 1:5) == Set([1, 2])

        sxcx = SimplicialComplex()
        push!(sxcx, (1, 3, 4, 5))
        push!(sxcx, (2, 3, 4, 5))
        push!(sxcx, (1, 3, 6))

        @test maximalcofaces(sxcx, (3, 4, 5)) == Set([1, 2])
        @test maximalcofaces(sxcx, (1, 3)) == Set([1, 3])
        @test maximalcofaces(sxcx, (1, 2, 3)) == Set()

        @test maximalfaces(sxcx, (1, 3, 4, 5, 6)) == Set([1, 3])
        @test maximalfaces(sxcx, (1, 2, 3, 4, 5)) == Set([1, 2])
        @test maximalfaces(sxcx, (1, 3, 5, 6)) == Set([3])
        @test maximalfaces(sxcx, (1, 2, 3)) == Set()
    end

    @testset "Insertion" begin
        sxcx = SimplicialComplex()
        push!(sxcx, (1, 3, 4, 5))
        push!(sxcx, (2, 3, 4, 5))
        push!(sxcx, (1, 3, 6))

        @test length(alllabels(sxcx)) == 3

        push!(sxcx, (1, 3, 4, 5))
        @test length(alllabels(sxcx)) == 3

        push!(sxcx, (1, 3, 4, 5, 7))
        @test length(alllabels(sxcx)) == 3

        push!(sxcx, (1, 3, 4, 5, 6, 7))
        @test length(alllabels(sxcx)) == 2

        push!(sxcx, (1, 2, 3, 4, 5, 6, 7))
        @test length(alllabels(sxcx)) == 1
    end

    @testset "in" begin
        sxcx = SimplicialComplex()
        push!(sxcx, (1, 3, 4, 5))
        push!(sxcx, (2, 3, 4, 5))
        push!(sxcx, (1, 3, 6))

        @test (1, 3, 4, 5) in sxcx
        @test (2, 3, 4, 5) in sxcx
        @test (1, 3, 6) in sxcx

        @test (1, 3) in sxcx
        @test (1, 6) in sxcx
        @test (2, 3, 5) in sxcx

        @test !((2, 3, 6) in sxcx)
        @test !((1, 7) in sxcx)

        sxcx = SimplicialComplex()
        push!(sxcx, (1, 2))
        push!(sxcx, (2, 3))
        push!(sxcx, (1, 3))
        @test !((1, 2, 3) in sxcx)

        @test_throws ArgumentError (2, 1) in sxcx
        @test_throws ArgumentError (7, 2) in sxcx
    end

    @testset "Removal" begin
        @testset "deletemaximals!" begin
            sxcx = SimplicialComplex()
            deletemaximals! = SimplicialComplexes.deletemaximals!
            push!(sxcx, (1, 2, 3))
            push!(sxcx, (2, 3, 4))
            push!(sxcx, (1, 2, 5))
            deletemaximals!(sxcx, Set([1, 2, 3]))

            @test isempty(sxcx)

            push!(sxcx, (1, 2, 3))
            push!(sxcx, (2, 3, 4))
            push!(sxcx, (1, 2, 5))

            deletemaximals!(sxcx, Set([1, 3]))
            @test length(SimplicialComplexes.maximalsimplices(sxcx)) == 1
        end
    end
end

@testset "Structure integrity" begin
    @testset "Labels" begin
        # Do stuff to a complex and make sure lbl_vals ∩ lbl_holes == {}
    end
end
