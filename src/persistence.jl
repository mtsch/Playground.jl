using LightGraphs
using SimpleWeightedGraphs

function getindexmatrix(dists)
    M = fill(0, size(dists))
    for (i, (u, v)) in enumerate(getsortededges(dists))
        M[u, v] = i
        M[v, u] = i
    end
    M
end

struct PersistenceState{T}
    # Precomputed properties.
    radius      ::T
    dists       ::Matrix{T}
    parents     ::Matrix{Int}
    indexmatrix ::Matrix{Int}

    # Internal state
    reduced    ::Vector{BitSet}
    generators ::Vector{Vector{Int}}

    # Buffers
    σ     ::BitSet
    cycle ::Vector{Int}
end

function PersistenceState(gc::GeodesicComplex)
    floyd = floyd_warshall_shortest_paths(gc)
    dists = floyd.dists
    parents = floyd.parents

    # todo: split connected components
    indexmatrix = getindexmatrix(dists)

    # buffers etc.
    reduced = fill(BitSet(), binomial(nv(gc), 2))
    generators = Vector{Int}[]

    PersistenceState(gc.radius, dists, parents, indexmatrix, reduced, generators, BitSet(), Int[])
end

function reduce!(st::PersistenceState)
end

# TODO
triangles(st::PersistenceState) = equilaterals(st.dists, st.radius)

# TODO: abstractgraph
function persistence(gc::GeodesicComplex; showprogress = false)
    showprogress && println("Calculating intrinsic persistence...")

    st = PersistenceState(gc)

    for tri in triangles(st)
        movecycle!(σ, cycle, parents, indexmatrix, tri)
        allunique(cycle) || continue

        showprogress && print(tri, " ↦ ", cycle)
        low = reduce!(σ, reduced)
        if low ≠ 0
            reduced[low] = copy(σ)
            if true || length(cycle) > 3
                push!(generators, copy(cycle))
            end
            showprogress && println(" survived.")
        else
            showprogress && println(" is kill.")
        end
        # check if minimal
    end
    generators
end

function reduce!(col, reduced)
    #println("Reduce $col")
    low = getlow(col)
    while low ≠ 0 && !isempty(reduced[low])
        #println("   ", col, " + ", reduced[low])
        col = sumcols!(col, reduced[low])
        low = getlow(col)
        #if length(col) == 3
        #    return 0
        #end
    end
    low
end

sumcols!(σ1::BitSet, σ2::BitSet) = setdiff!(σ1, σ2)
getlow(σ::BitSet) = length(σ) > 0 ? last(σ) : 0

# TODO: get points and edges
# TODO: islocallyminimal
function movecycle!(σ, cycle, parents, indexmatrix, (i, j, k))
    empty!(σ)
    empty!(cycle)

    addpath!(σ, cycle, parents, indexmatrix, i, j)
    addpath!(σ, cycle, parents, indexmatrix, j, k)
    addpath!(σ, cycle, parents, indexmatrix, k, i)
end

function addpath!(σ, cycle, parents, indexmatrix, v, u)
    while v ≠ u
        v′ = parents[u, v]
        push!(σ, indexmatrix[v, v′])
        push!(cycle, v)
        v = v′
    end
end
