sal = SAL()

push!(sal, (1,3,4,5))
push!(sal, (2,3,4,5))
push!(sal, (1,3,6))
println(values(sal.nodelist[1]))

println(collect(SimplexTrees.maximalsimplices(sal)))

@testset "SAL: DictSet" begin
    ds = SimplexTrees.DictSet{Int}()

    # Insertions.
    insert!(ds, 1, 1)
    insert!(ds, 1, 2)
    insert!(ds, 1, 3)
    @test ds[1] == IntSet([1, 2, 3])
    @test ds[2] == IntSet()
    @test ds[3] == IntSet()

    insert!(ds, 2, 2)
    @test ds[1] == IntSet([1, 2, 3])
    @test ds[2] == IntSet(2)
    @test ds[3] == IntSet()

    # Deletions, isempty
    delete!(ds, 1)
    @test ds[1] == IntSet()
    @test ds[2] == IntSet(2)

    insert!(ds, 1, 1)
    insert!(ds, 1, 2)
    insert!(ds, 1, 3)
    delete!(ds, 1, 2)
    @test ds[1] == IntSet([1, 3])
    @test ds[2] == IntSet(2)
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
    delete!(ds, IntSet([2, 3]))
    @test ds[1] == IntSet(1)
    @test ds[2] == IntSet()
    @test ds[3] == IntSet()

    delete!(ds, IntSet([1, 2, 3]))
    @test isempty(ds)
end

@testset "isface" begin
    @test isface((1, 2, 3), (1, 2, 3))
    @test isface([3, 6, 7], 1:10)
    @test isface([70, 80], 1:100)
    @test !isface([70, 80], [1])
    @test !isface([1, 6, 7], 1:6)
    @test !isface([1, 6, 7], [1, 2, 3, 4, 5, 7, 8, 9])
    @test !isface([1, 6, 7], [1, 7])
end

@testset "maximalsimplices" begin
    sal = SAL()
    push!(sal, (1, 2, 3, 4))
    push!(sal, (3, 4, 5))
    maxsxs = collect(values(SimplexTrees.maximalsimplices(sal)))
    @test [3, 4, 5] in maxsxs
    @test [1, 2, 3, 4] in maxsxs
    @test !([1, 2, 3] in maxsxs)
end

@testset "Labels" begin
    # Make sure lbl_vals ∩ lbl_holes == {}
end

@testset "SAL: in" begin
    sal = SAL()
    push!(sal, (1,3,4,5))
    push!(sal, (2,3,4,5))
    push!(sal, (1,3,6))

    @test (1, 3, 4, 5) in sal
    @test (2, 3, 4, 5) in sal
    @test (1, 3, 6) in sal

    @test (1, 3) in sal
    @test (1, 6) in sal
    @test (2, 3, 5) in sal

    @test !((2, 3, 6) in sal)
    @test !((1, 7) in sal)

    @test_throws ArgumentError (2, 1) in sal
    @test_throws ArgumentError (7, 2) in sal

    sal = SAL()
    push!(sal, (1, 2))
    push!(sal, (2, 3))
    push!(sal, (1, 3))
    @test !((1, 2, 3) in sal)
end

@testset "SAL: Insertion" begin
end

@testset "SAL: Removal" begin
    @testset "removemaximals" begin
        sal = SAL()
        push!(sal, (1, 2, 3))
        push!(sal, (2, 3, 4))
        push!(sal, (1, 2, 5))
        SimplexTrees.deletemaximals!(sal, IntSet([1, 2, 3]))

        @test isempty(sal)

        push!(sal, (1, 2, 3))
        push!(sal, (2, 3, 4))
        push!(sal, (1, 2, 5))

        SimplexTrees.deletemaximals!(sal, IntSet([1, 3]))
        @test length(SimplexTrees.maximalsimplices(sal)) == 1
    end
end
