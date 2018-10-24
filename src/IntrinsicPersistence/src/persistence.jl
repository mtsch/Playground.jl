# In M, a n×n distance matrix, find s such that:
# ∀i ∈ 1:n, ∃j ∈ 1:n, i≠j; M[i, j] ≤ s.
function sampledensity(M::AbstractMatrix{T}) where T
    s = typemin(T)
    for i in 1:size(M, 1)
        nearest = typemax(T)
        for j in 1:size(M, 1)
            i == j && continue
            dist = M[j, i]
            nearest = dist < nearest ? dist : nearest
        end
        s = s > nearest ? s : nearest
    end
    s
end

sampledensity(g::AbstractGraph) = maximum(weights(g))

function edgetoorder(G::AbstractGraph)
    w = adjacency_matrix(G)
    for (i, e) in enumerate(edges(G))
        w[src(e), dst(e)] = i
        w[dst(e), src(e)] = i
    end
    Symmetric(w)
end
function edgetoorder(G::AbstractSimpleWeightedGraph)
    es = collect(edges(G))
    sort!(es, by = weight)
    I = Vector{Int}(length(es))
    J = similar(I)
    i = 1
    for e in es
        I[i] = src(e)
        J[i] = dst(e)
        i += 1
    end
    n = nv(G)
    Symmetric(sparse(I, J, 1:length(es), n, n), :U)
end

function spannergraph(pts::AbstractMatrix{T},
                      metric = Euclidean(), r = nothing
                      ) ::SimpleWeightedGraph{Int, T} where T
    s    = r == nothing ? 2 * sampledensity(pairwise(metric, pts)) : T(r)
    n, m = size(pts)
    _pts = unique(reinterpret(SVector{n, T}, pts, (m,)))
    spannergraph(_pts, metric, s)
end
function spannergraph(pts::AbstractVector{SVector{D, T}},
                      metric = Euclidean(), r = nothing
                      ) ::SimpleWeightedGraph{Int, T} where {D, T}
    n = length(pts)
    s = r == nothing ?
        2 * T(sampledensity(pairwise(metric, reinterpret(T, pts, (D, n))))) :
        T(r)
    btree  = BallTree(pts, metric, reorder = true)
    neighs = inrange(btree, pts, s, true)
    is = Int[]; js = Int[]; ds = T[]
    for (i, ns) in enumerate(neighs)
        for j in ns
            i < j || continue
            append!(is, [i, j])
            append!(js, [j, i])
            d = evaluate(metric, pts[i], pts[j])
            append!(ds, [d, d])
        end
    end
    n = length(pts)
    SimpleWeightedGraph(sparse(is, js, ds, n, n))
end

function persistence(pts::AbstractVecOrMat;
                     show_trace = false,
                     use_hash   = false,
                     skip_small = true,
                     resolution = nothing,
                     metric     = Euclidean())

    graph = spannergraph(pts, metric, resolution)
    res   = persistence(graph,
                        eltype(weights(graph)),
                        show_trace = show_trace,
                        skip_small = skip_small,
                        use_hash   = use_hash)
    add_data(res, pts)
end

function persistence(graph::AbstractGraph,
                     T::DataType = eltype(weights(graph));
                     skip_small ::Bool = true,
                     show_trace ::Bool = false,
                     use_hash   ::Bool = false) # <- T mora bit eltype(weights(G))
                                                        #    preveri če mora bit parameter

    is_directed(graph) && throw(ArgumentError("Directed graphs not supported!"))

    show_trace && print("# ===================== #\n",
                        "| Intrinsic persistence |\n",
                        "# ===================== #\n\n",
                        "Preprocessing...")

    s = sampledensity(graph)
    n = nv(graph)

    # Constants.
    order     = edgetoorder(graph)
    dijkstra  = parallel_multisource_dijkstra_shortest_paths(graph)
    distmx    = dijkstra.dists
    parents   = dijkstra.parents
    nrows     = ne(graph)
    triangles = equilaterals(distmx, tol = s)

    # Result.
    generators = IntrinsicGenerator{T}[]

    # State.
    reductionmx      = fill(Int[], nrows)
    st               = start(triangles)
    ngenerators_left = nrows - n + length(connected_components(graph))

    # Debug variables.
    r = zero(T)
    currindex = 0
    nskipped  = 0
    nreduced  = 0

    # Buffers.
    column = Int[]; sizehint!(column, n)
    buff   = Int[]; sizehint!(buff, n)
    paths  = Int[], Int[], Int[]

    if use_hash
        hashes = Set{UInt}()
    end

    if show_trace
        @time ntriangles = iterlength(triangles)
        print("\r")
        println("* Sample density:       $(round(s, 5))")
        println("* Spanner density:      $(density(graph))")
        println("* Connected components: $(length(connected_components(graph)))")
        println("* Number of triangles:  $(ntriangles)")
        println("* Number of generators: $(ngenerators_left)")
        println()
    end

    maxcollen = 0
    while ngenerators_left > 0 && !done(triangles, st)
        currindex += 1
        Δ, st = next(triangles, st)
        i, j, k = Δ.vertices
        r = Δ.r

        if show_trace && currindex % (ntriangles ÷ 100) == 0
            showprogress(currindex, r, ntriangles, nreduced,
                         length(generators), ngenerators_left, nrows)
        end

        empty!(column); empty!(paths[1]); empty!(paths[2]); empty!(paths[3])
        collectpath!(column, paths[1], i, j, parents, order)
        collectpath!(column, paths[2], j, k, parents, order)
        collectpath!(column, paths[3], k, i, parents, order)
        sort!(column)

        # Skipping.
        #r < 2s && length(column) > 3           && continue
        use_hash && !add_hash!(hashes, column) && continue
        !alluniquesorted(column)               && continue
        !islocallyminimal(paths, distmx, s)    && continue

        # Matrix reduction.
        low = getlow(column)
        stop = false
        while !stop && low ≠ 0 && !isempty(reductionmx[low])
            column, buff = sumcols!(buff, column, reductionmx[low])
            low = getlow(column)
            maxcollen = max(length(column), maxcollen)

            if use_hash && !add_hash!(hashes, column)
                stop = true
            end
        end
        stop && continue

        # Output.
        if low ≠ 0
            reductionmx[low] = column
            column = Int[]

            if (!skip_small && sum(length.(paths)) > 3) || r > s
                path = Int[]
                for i in 1:3
                    append!(path, paths[i])
                end
                lengths = similar(path, T)
                for p_i in eachindex(path)
                    p_j = p_i < length(path) ? p_i + 1 : 1
                    lengths[p_i] = distmx[path[p_i], path[p_j]]
                end

                push!(generators, IntrinsicGenerator(Δ, path, lengths))
            end
            ngenerators_left -= 1
            nskipped += 1
        else
            nreduced += 1
        end
    end

    if show_trace
        showprogress(currindex, r, ntriangles, nreduced,
                     length(generators), ngenerators_left, nrows)
        if use_hash
            println("no. hashes: $(length(hashes))")
        end
        println()
        println()
        @show maxcollen
    end
    ngenerators_left > 0 && warn("Not all generators found!")

    IntrinsicPersistenceResult(generators, graph)
end

function collectpath!(∂_j, side, v1, v2, parents, order)
    # sem not gre lahko čekiranje če je prejšnji trikotnik vsebovan v poti
    index = v1
    while parents[v2, index] ≠ 0
        push!(side, index)
        nxtindex = parents[v2, index]
        push!(∂_j, order[index, nxtindex])
        index = nxtindex
    end
end

function add_hash!(hashes, column)
    h = hash(column)
    if h in hashes
        false
    else
        push!(hashes, h)
        true
    end
end

function alluniquesorted(vec)
    for i in 1:length(vec)-1
        vec[i] == vec[i+1] && return false
    end
    true
end

# Check for each pair of points on opposing sides of the triangle whether the
# path through triangle is the shortest path.
function islocallyminimal(sides, distmx, s)
    @inbounds for (i,j,k) in ((1,2,3), (2,3,1), (3,1,2))
        a = sides[i]; b = sides[j]; c = sides[k]
        common = b[1]
        opp1   = a[1]
        opp2   = c[1]

        for i in 2:length(a), j in 2:length(b)
            direct = distmx[a[i], b[j]]
            around = min(distmx[a[i], common] + distmx[common, b[j]],
                         distmx[a[i], opp1] + distmx[opp1, opp2] + distmx[opp2, b[j]])
            if direct + s < around
                return false
            end
        end
    end
    true
end

function getlow(col::AbstractVector{T}) where T
    isempty(col) ? zero(T) : last(col)
end

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

function iterlength(itr)
    len = 0
    for _ in itr
        len += 1
    end
    len
end

function showprogress(j, r, ntriangles, nreduced, ngenerators, ngenerators_left, nrows)
    print("\rTriangles visited:",
          lpad(round(Int, 100 * j/ntriangles), 3), "%, ℓ ≈ ",
          rpad(round(3r, 5), 7, 0),
          ", reduced:",
          lpad(round(Int, 100 * nreduced/j), 3), "%",
          ". Generators found:",
          lpad(ngenerators, 2),
          ", left: ",
          rpad(ngenerators_left, floor(Int, log10(nrows) + 1)))
    flush(STDOUT::IO)
end
