"""
    alluniquesorted(vec)

Check if all values of sorted vector `vec` are unique.
"""
function alluniquesorted(vec)
    for i in 1:length(vec)-1
        vec[i] == vec[i+1] && return false
    end
    true
end

"""
    sumcols!(buff, col1, col2)

Calculate the symmetric difference of two sorted integer vectors. This is
used to sum columns in the reduction matrix (mod 2).

Writes the result to `buff` and also returns an emptied `col1`, which can be
reused as `buff` in the next iteration.
"""
@inline function sumcols!(buff, col1, col2)
    i = j = 1
    n = length(col1)
    m = length(col2)
    @inbounds while i ≤ n && j ≤ m
        if col1[i] < col2[j]
            push!(buff, col1[i])
            i += 1
        elseif col1[i] > col2[j]
            push!(buff, col2[j])
            j += 1
        else
            i += 1
            j += 1
        end
    end
    @inbounds while j ≤ m
        push!(buff, col2[j])
        j += 1
    end
    @inbounds while i ≤ n
        push!(buff, col1[i])
        i += 1
    end
    empty!(col1)
    buff, col1
end

"""
    getlow(vec)

Get the "low one" in a column of the reduction matrix.
"""
function getlow(col::AbstractVector{T}) where T
    isempty(col) ? zero(T) : last(col)
end


"""
    allpaths(G::SimpleWeightedGraph)

Get all shortest paths in graph.
"""
function allpaths(G::SimpleWeightedGraph{I, T}) where {I, T}
    paths = Vector{DijkstraState{T, I}}(nv(G))
    for v in vertices(G)
        paths[v] = dijkstra_shortest_paths(G, v)
    end
    paths
end

"""
    iterlength(itr)

Get the length of an iterator that does not have `length` defined.
"""
function iterlength(itr)
    len = 0
    for _ in itr
        len += 1
    end
    len
end

"""
    addpath!(res_e, res_p, v1, v2, paths, edgeord)

Find the path from `v1` to `v2` in `paths`, convert to edge indices from
`edgeord` and append the result to `res`.
"""
function addpath!(res_e, res_p, v1, v2, paths, edgeord)
    @inbounds begin
        path = enumerate_paths(paths[v1], v2)
        for i in 1:length(path)-1
            push!(res_e, edgeord[path[i], path[i + 1]])
            push!(res_p, path[i])
        end
    end
    res_e, res_p
end

# Rewrite
@inline function islocallyminimal(sides, distmat, s)
    s12, s23, s31 = sides
    l12, l23, l31 = length.(sides)
    for i in 2:l12, j in 2:l23
        if distmat[s12[i], s23[j]] < distmat[s12[i], s23[1]] + distmat[s23[1], s23[j]] - 2s
            return false
        end
    end
    for i in 2:l23, j in 2:l31
        if distmat[s23[i], s31[j]] < distmat[s23[i], s31[1]] + distmat[s31[1], s31[j]] - 2s
            return false
        end
    end
    for i in 2:l31, j in 2:l12
        if distmat[s31[i], s12[j]] < distmat[s31[i], s12[1]] + distmat[s12[1], s12[j]] - 2s
            return false
        end
    end
    true
end

function showprogress(j, r, ntriangles, nskipped, ngenerators, ngenerators_left, nrows)
    print("\rTriangles visited:",
          lpad(round(Int, 100 * j/ntriangles), 3), "%, l ≈ ",
          rpad(round(3r, 5), 7, 0),
          ", skipped:",
          lpad(round(Int, 100 * nskipped/j), 3), "%",
          ". Generators found:",
          lpad(ngenerators, 2),
          ", left: ",
          rpad(ngenerators_left, floor(Int, log10(nrows) + 1)))
    flush(STDOUT::IO)
end

"""
    persistence(points; printprogress, resolution, subset)

Calculate 1-dimensional intrinsic persistence of `points`.
"""
function persistence(pts::AbstractMatrix{T}; kwargs...) where T
    n, m = size(pts)
    _pts = reinterpret(SVector{n, T}, pts, (m,))
    persistence(_pts; kwargs...)
end
function persistence(pts::AbstractVector{SVector{D, T}};
                     show_trace = false,
                     resolution = nothing) where {D, T}

    # Prepare default arguments.
    if resolution == nothing
        s = sampledensity(pairwise(Euclidean(),
                                   reinterpret(T, pts, (D, length(pts)))))
    else
        s = resolution
    end

    # Main algo.
    show_trace && print("# ===================== #\n",
                        "| Intrinsic persistence |\n",
                        "# ===================== #\n\n",
                        "Preprocessing...")

    # Intrinsic metric distance matrix and adjacency graph.
    metric  = IntrinsicMetric(pts, 2s)
    distmat = pairwise(metric)
    graph   = adjgraph(metric)

    # Number of rows in reduction matrix.
    nrows = ne(graph)
    # Edge ordering matrix.
    edgeord = edgetoorder(graph)
    # Shortest paths from every node.
    paths = allpaths(graph)
    # Triangle iterator.
    triangles = equilaterals(distmat, tol = 2s, lower = 0)
    # Reduction matrix columns indexed by getlow.
    L = fill(Int[], nrows)

    # Result storage.
    generators = IntrinsicGenerator{T}[]
    # Number of generators left to discover.
    ngenerators_left = nrows - length(pts) + length(connected_components(graph))
    nskipped = 0

    if show_trace
        ntriangles = iterlength(triangles)
        print("\r")
        println("* Sample density:       $(round(s, 5))")
        println("* Number of triangles:  $ntriangles")
        println("* Number of generators: $(ngenerators_left)")
        println()
    end

    # j-th column and buffer are swapped back and forth to avoid allocation.
    ∂_j   = Int[]
    buff  = Int[]
    # Sides of the triangle as paths.
    sides = (Int[], Int[], Int[])
    currindex = 0
    r = zero(T)
    for Δ in triangles
        currindex += 1
        ngenerators_left > 0 || break

        v1, v2, v3 = Δ.vertices
        r = Δ.r

        if show_trace && currindex % (ntriangles ÷ 100) == 0
            showprogress(currindex, r, ntriangles, nskipped,
                         length(generators), ngenerators_left, nrows)
        end

        # Build column and path.
        empty!(∂_j)
        empty!(sides[1]); empty!(sides[2]); empty!(sides[3])
        addpath!(∂_j, sides[1], v1, v2, paths, edgeord)
        addpath!(∂_j, sides[2], v2, v3, paths, edgeord)
        addpath!(∂_j, sides[3], v3, v1, paths, edgeord)
        sort!(∂_j)

        # Skip
        #!alluniquesorted(∂_j) && continue

        if true
        if (r ≤ 2s && length(∂_j) > 3) ||
           !alluniquesorted(∂_j) ||
           !islocallyminimal(sides, distmat, s)
            nskipped += 1
            continue
        end
        end

        # Matrix reduction.
        low = getlow(∂_j)
        while low ≠ 0 && !isempty(L[low])
            ∂_j, buff = sumcols!(buff, ∂_j, L[low])
            low = getlow(∂_j)
        end

        # Output.
        if low ≠ 0
            L[low] = ∂_j
            ∂_j = Int[]

            if r > 2s
                path = Int[]
                for i in 1:3
                    append!(path, sides[i])
                end
                lengths = similar(path, T)
                for p_i in eachindex(path)
                    p_j = p_i < length(path) ? p_i + 1 : 1
                    lengths[p_i] = distmat[path[p_i], path[p_j]]
                end

                push!(generators, IntrinsicGenerator(Δ, path, lengths))
            end

            ngenerators_left -= 1
            nskipped += 1
        end
    end
    if show_trace
        showprogress(currindex, r, ntriangles, nskipped,
                     length(generators), ngenerators_left, nrows)
        println()
        println()
    end

    IntrinsicPersistenceResult(generators, pts, collect(1:length(pts)))
end
