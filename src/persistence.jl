using LightGraphs
using SimpleWeightedGraphs

getradius(g::AbstractGraph) = maximum(weights(g))
getradius(g::GeodesicComplex) = g.radius

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

function PersistenceState(gc::AbstractGraph)
    floyd = floyd_warshall_shortest_paths(gc)
    dists = floyd.dists
    parents = floyd.parents

    # todo: split connected components
    indexmatrix = getindexmatrix(dists)

    # buffers etc.
    reduced = fill(BitSet(), binomial(nv(gc), 2))
    generators = Vector{Int}[]

    PersistenceState(getradius(gc), dists, parents, indexmatrix, reduced, generators, BitSet(), Int[])
end

Base.show(io::IO, st::PersistenceState) = print(io, "PersistenceState")

#sumcols!(σ1::BitSet, σ2::BitSet) = symdiff!(σ1, σ2)
getlow(σ::BitSet) = length(σ) > 0 ? last(σ) : 0

function reduce!(st::PersistenceState)
    low = getlow(st.σ)
    while low ≠ 0 && !isempty(st.reduced[low])
        symdiff!(st.σ, st.reduced[low])
        low = getlow(st.σ)
    end
    low
end

# TODO
eqtriangles(st::PersistenceState) = equilaterals(st.dists, 2st.radius)

# TODO: get points and edges
# TODO: islocallyminimal
function movecycle!(st, (i, j, k))
    empty!(st.σ)
    empty!(st.cycle)

    addpath!(st, i, j)
    addpath!(st, j, k)
    addpath!(st, k, i)
end

function addpath!(st, v, u)
    while v ≠ u
        v′ = st.parents[u, v]
        push!(st.σ, st.indexmatrix[v, v′])
        push!(st.cycle, v)
        v = v′
    end
end

# TODO: abstractgraph
function persistence(gc::AbstractGraph; showprogress = false)
    showprogress && println("Calculating intrinsic persistence...")

    st = PersistenceState(gc)

    for (i, tri) in enumerate(eqtriangles(st))
        movecycle!(st, tri)
        length(st.σ) == length(st.cycle) || continue

        showprogress && print(tri, " ↦ ", st.cycle, ", ", collect(st.σ))
        low = reduce!(st)
        if low ≠ 0
            st.reduced[low] = copy(st.σ)
            if true || length(st.cycle) > 3
                push!(st.generators, copy(st.cycle))
            end
            showprogress && println(" survived as ", collect(st.σ), ".")
        else
            showprogress && println(" is kill.")
        end
    end
    st
end
