function test_structure(sxt::SXTree, pairs...)
    function count_children(sxt, seq)
        curr = sxt[seq[1]]
        for i in seq[2:end]
            curr = curr[i]
        end
        SimplexTrees.num_children(curr)
    end

    all(count_children(sxt, seq) == n for (seq, n) in pairs)
end

function iter_count(f, n)
    i = 0
    for n in f(n)
        i += 1
    end
    i
end

function iter_collect(f, n)
    res = Vector{}[]
    for n in f(n)
        push!(res, vec(n))
    end
    res
end

@testset "Simplex insertion" begin
    @testset "Single simplex" begin

        # Check the structure of a simplex.
        # Build a 4-simplex and test the structure.
        structure4(sxt) = test_structure(sxt,
                                         [1] => 3,
                                         [1, 2] => 2,
                                         [1, 2, 3] => 1,
                                         [1, 2, 3, 4] => 0,
                                         [1, 2, 4] => 0,
                                         [1, 3] => 1,
                                         [1, 3, 4] => 0,

                                         [2] => 2,
                                         [2, 3] => 1,
                                         [2, 3, 4] => 0,
                                         [2, 4] => 0,

                                         [3] => 1,
                                         [3, 4] => 0,

                                         [4] => 0)
        sxt = SXTree(4)
        insert_simplex!(sxt, collect(1:4))
        @test structure4(sxt)

        # Build a full 6-simplex and roughly test the structure.
        structure6(sxt) = test_structure(sxt,
                                         Pair.([1:i for i in 1:6], 5:-1:0)...,
                                         Pair.([2:i for i in 2:6], 4:-1:0)...,
                                         Pair.([3:i for i in 3:6], 3:-1:0)...,
                                         Pair.([4:i for i in 4:6], 2:-1:0)...,
                                         Pair.([4:i for i in 5:6], 1:-1:0)...,
                                         [6] => 0)

        sxt = SXTree(6)
        insert_simplex!(sxt, collect(1:6))
        @test structure6(sxt)

    end

    @testset "Case from article" begin

        # Check the structure of a complex.
        # See: https://hal.archives-ouvertes.fr/file/index/docid/707901/filename/RR-7993.pdf
        # page 5, first connected component for the image.
        article_structure(sxt) = test_structure(sxt,
                                                [1] => 2,
                                                [1, 2] => 1,
                                                [1, 2, 3] => 0,
                                                [1, 3] => 0,

                                                [2] => 3,
                                                [2, 3] => 2,
                                                [2, 3, 4] => 1,
                                                [2, 3, 4, 5] => 0,
                                                [2, 4] => 1,
                                                [2, 4, 5] => 0,
                                                [2, 5] => 0,

                                                [3] => 2,
                                                [3, 4] => 1,
                                                [3, 4, 5] => 0,
                                                [3, 5] => 0,

                                                [4] => 1,
                                                [4, 5] => 0,

                                                [5] => 0)

        # Insert some simplices and make sure all of them were inserted properly.
        sxt = SXTree(5)
        # Sanity check.
        @test n_vertices(sxt) == 5

        insert_simplex!(sxt, [1, 3, 2])
        insert_simplex!(sxt, [2, 3, 4, 5])

        @test article_structure(sxt)
        @test dim(sxt) == 3

        # Insert in a different order
        sxt = SXTree(5)
        insert_simplex!(sxt, [1, 2])
        insert_simplex!(sxt, [5, 4])
        insert_simplex!(sxt, [3, 2])
        insert_simplex!(sxt, [3, 5])
        insert_simplex!(sxt, [1, 3])
        insert_simplex!(sxt, [2, 1, 3])
        insert_simplex!(sxt, [4, 5, 3, 2])

        @test article_structure(sxt)
        @test dim(sxt) == 3

        # Add some simplices that were already added and check again.
        insert_simplex!(sxt, [3, 4])
        insert_simplex!(sxt, [3, 4])
        insert_simplex!(sxt, [3, 4])
        insert_simplex!(sxt, [3, 4])

        @test article_structure(sxt)

        # Invalid arguments.
        @test_throws ArgumentError insert_simplex!(sxt, [1, 5, 6])
        @test_throws ArgumentError insert_simplex!(sxt, [1, 5, 2, 3, 1])

        # Nothing broke when attempting to call with invalid arguments.
        @test article_structure(sxt)
        @test dim(sxt) == 3
    end
end


@testset "Iterators" begin

    @testset "Cousins" begin
        cousin_count(sxt) = iter_count(cousins, sxt)
        # One element list.
        sxt = SXTree(2)
        insert_simplex!(sxt, [1, 2])

        @test first(cousins(sxt[1])) == sxt[1]
        @test cousin_count(sxt[1][2]) == 1

        # Build the test complex.
        sxt = SXTree(5)
        insert_simplex!(sxt, [1, 3, 2])
        insert_simplex!(sxt, [2, 3, 4, 5])

        @test first(cousins(sxt[1])) == sxt[1]

        @test cousin_count(sxt[2, 5]) == cousin_count(sxt[4][5]) == 3
        @test cousin_count(sxt[2, 3, 5]) == cousin_count(sxt[3][4][5]) == 3
    end

    @testset "Upwards" begin
        up_count(sxt) = iter_count(ancestors, sxt)

        sxt = SXTree(6)
        insert_simplex!(sxt, collect(1:6))

        @test up_count(sxt[6]) == 1
        @test up_count(sxt[1, 2, 3, 4, 5, 6]) == 6

        @test label.(collect(ancestors(sxt[4, 5, 6]))) == [6, 5, 4]

        @test length(ancestors(sxt[1, 2, 3, 4, 5, 6])) == 6
    end

    @testset "Vector representation" begin
        sxt = SXTree(6)
        insert_simplex!(sxt, collect(1:6))

        @test vec(sxt[1, 2]) == [1, 2]
        @test vec(sxt[2, 5, 6]) == [2, 5, 6]
        @test vec(sxt[1, 2, 3, 4, 5, 6]) == [1, 2, 3, 4, 5, 6]
    end

    @testset "BF iteration" begin
        sxt = SXTree(3)

        succ_bf_count(sxt) = iter_count(successors, sxt)
        succ_bf_collect(sxt) = iter_collect(successors, sxt)

        @test succ_bf_count(sxt) == 3
        @test succ_bf_count(sxt[3]) == 0
        insert_simplex!(sxt, collect(1:3))
        @test succ_bf_count(sxt) == 7
        @test succ_bf_collect(sxt) == [[1], [2], [3], [1,2], [1,3], [2,3], [1,2,3]]
        @test succ_bf_collect(sxt[1]) == [[1,2], [1,3], [1,2,3]]
    end

    @testset "DF iteration" begin
        sxt = SXTree(3)

        succ_df_count(sxt) = iter_count(x -> successors(x, DFIter), sxt)
        succ_df_collect(sxt) = iter_collect(x -> successors(x, DFIter), sxt)

        @test succ_df_count(sxt) == 3
        @test succ_df_count(sxt[3]) == 0
        insert_simplex!(sxt, collect(1:3))
        @test succ_df_count(sxt) == 7
        @test succ_df_collect(sxt) == [[1], [1,2], [1,2,3], [1,3], [2], [2,3], [3]]
        @test succ_df_collect(sxt[1]) == [[1,2], [1,2,3], [1,3]]
    end
end

@testset "Cofaces" begin
    sxt6 = SXTree(6)
    insert_simplex!(sxt6, collect(1:6))

    sxt4 = SXTree(4)
    insert_simplex!(sxt4, collect(1:4))

    @testset "is_face" begin
        @test is_face([1, 2, 3], sxt6[1, 2, 3, 4, 5, 6])
        @test is_face([4, 2, 6], sxt6[1, 2, 3, 4, 5, 6])
        @test !is_face([4, 2, 6], sxt6[1, 2, 3, 4, 5])

        @test_throws ArgumentError !is_face([1, 1], sxt6[1, 2, 3, 4, 5, 6])
    end

    @testset "cofaces" begin
        @test length(cofaces(sxt6, [1, 2, 3, 4, 5])) == 2
        @test length(cofaces(sxt6, [2, 3, 4, 5])) == 4
        @test length(cofaces(sxt6, [2, 4, 5])) == 8

        cfcs = vec.(cofaces(sxt4, [1, 2]))
        @test length(cfcs) == 4
        @test [1, 2] in cfcs &&
            [1, 2, 4] in cfcs &&
            [1, 2, 3] in cfcs &&
            [1, 2, 3, 4] in cfcs
    end
end

@testset "Simplex removal" begin
    @testset "rem_node!" begin
        sxt6 = SXTree(6)
        insert_simplex!(sxt6, collect(1:6))

        SimplexTrees.rem_node!(sxt6, sxt6[1,2,3,4,5,6])
        @test iter_count(successors, sxt6[1,2,3,4,5]) == 0

        sxt4 = SXTree(4)
        insert_simplex!(sxt4, collect(1:4))
        SimplexTrees.rem_node!(sxt4, sxt4[1, 2])
        @test iter_count(successors, sxt4[1]) == 3

        @test_throws ArgumentError SimplexTrees.rem_node!(sxt4, sxt4[1])
    end

    @testset "rem_simplex!" begin
        sxt6 = SXTree(6)
        insert_simplex!(sxt6, collect(1:6))

        rem_simplex!(sxt6, [1,2,3,4,5,6])
        @test dim(sxt6) == 4

        insert_simplex!(sxt6, collect(1:6))
        rem_simplex!(sxt6, sxt6[1,2,3,4,5,6])
        @test dim(sxt6) == 4

        rem_simplex!(sxt6, [1,2])
        @test dim(sxt6) == 4

        rem_simplex!(sxt6, [2,3])
        rem_simplex!(sxt6, [3,4])
        rem_simplex!(sxt6, [4,5])
        rem_simplex!(sxt6, [5,6])
        @test dim(sxt6) == 2
        @test n_simplices(sxt6) == 20
    end
end
