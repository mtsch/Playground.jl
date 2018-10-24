@testset "Utils" begin
    @testset "DictSet" begin
        ds = SimplicialComplexes.DictSet{Int, Int}()

        # Insertions.
        insert!(ds, 1, 1)
        insert!(ds, 1, 2)
        insert!(ds, 1, 3)
        @test ds[1] == Set{Int}([1, 2, 3])
        @test ds[2] == Set{Int}()
        @test ds[3] == Set{Int}()

        insert!(ds, 2, 2)
        @test ds[1] == Set{Int}([1, 2, 3])
        @test ds[2] == Set{Int}(2)
        @test ds[3] == Set{Int}()

        # Deletions, isempty
        delete!(ds, 1)
        @test ds[1] == Set{Int}()
        @test ds[2] == Set{Int}(2)

        insert!(ds, 1, 1)
        insert!(ds, 1, 2)
        insert!(ds, 1, 3)
        delete!(ds, 1, 2)
        @test ds[1] == Set{Int}([1, 3])
        @test ds[2] == Set{Int}(2)
        @test !isempty(ds)

        delete!(ds, 1, 1)
        delete!(ds, 1, 3)
        delete!(ds, 2)
        @test isempty(ds)

        # Deletions with sets
        insert!(ds, 1, 1)
        insert!(ds, 1, 2)
        insert!(ds, 1, 3)
        insert!(ds, 2, 3)
        insert!(ds, 3, 2)
        delete!(ds, Set{Int}([2, 3]))
        @test ds[1] == Set{Int}(1)
        @test ds[2] == Set{Int}()
        @test ds[3] == Set{Int}()

        delete!(ds, Set{Int}([1, 2, 3]))
        @test isempty(ds)
    end

    @testset "validatesx" begin
        validatesx = SimplicialComplexes.validatesx
        @test_throws ArgumentError validatesx((1, 1, 2))
        @test_throws ArgumentError validatesx((1, 3, 2))
        @test_throws ArgumentError validatesx([1, 5, 5])
        @test_throws ArgumentError validatesx([1, 5, 4])
        @test_throws ArgumentError validatesx([1.0, 2.0, 3.0])
        @test_throws ArgumentError validatesx(["foo", "bar", "baz"])
        @test_throws ArgumentError validatesx((1, 2, 3, '!'))
        @test_throws ArgumentError validatesx(Set([1, 2, 3]))
        @test validatesx(Int[]) == []
        @test validatesx(()) == ()
        @test validatesx(1:10) == 1:10
        @test validatesx((1, 2, 3)) == (1, 2, 3)
        @test validatesx(Int32[1, 200, 3000]) == [1, 200, 3000]
    end

    @testset "isface" begin
        isface = SimplicialComplexes.isface
        @test isface((1, 2, 3), (1, 2, 3))
        @test isface([3, 6, 7], 1:10)
        @test isface([70, 80], 1:100)
        @test !isface([70, 80], [1])
        @test !isface([1, 6, 7], 1:6)
        @test !isface([1, 6, 7], [1, 2, 3, 4, 5, 7, 8, 9])
        @test !isface([1, 6, 7], [1, 7])
    end
end
