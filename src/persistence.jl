using LightGraphs
using SimpleWeightedGraphs
using StaticArrays
using Distances
using NearestNeighbors

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
    dists       ::Matrix{T}
    parents     ::Matrix{Int}
    indexmatrix ::Matrix{Int}
    triangles   ::EquilateralIterator{T}
    # Internal state
    reduced     ::Vector{BitSet}
    generators  ::Vector{Vector{Int}}
    diameters   ::Vector{T}
    # Buffers
    σ           ::BitSet
    cycle       ::Vector{Int}
end

coverdiameter(g::AbstractGraph) = maximum(weights(g))
coverdiameter(g::GeodesicComplex) = 2g.radius

function PersistenceState(gc::AbstractGraph)
    d = coverdiameter(gc)
    # todo: split connected components
    floyd = floyd_warshall_shortest_paths(gc)
    dists = floyd.dists
    parents = floyd.parents
    indexmatrix = getindexmatrix(dists)
    triangles = equilaterals(dists, d)

    reduced = fill(BitSet(), binomial(nv(gc), 2))

    PersistenceState(dists, parents, indexmatrix, triangles,
                     reduced, Vector{Int}[], eltype(dists)[], BitSet(), Int[])
end

Base.show(io::IO, st::PersistenceState) = print(io, "PersistenceState")

"""
    low(σ::BitSet)

Get the largest index in `σ`.
"""
low(σ::BitSet) = length(σ) > 0 ? last(σ) : 0

"""
    reduce!(st::PersistenceState)

Use the standard reduction algorithm to reduce `st.σ`.
"""
function reduce!(st::PersistenceState)
    l = low(st.σ)
    while l ≠ 0 && !isempty(st.reduced[l])
        symdiff!(st.σ, st.reduced[l])
        l = low(st.σ)
    end
    l
end

"""
    movecycle!(st::PersistenceState, Δ)

Move cycle spanned by `Δ` into `st.cycle` and `st.σ`.
"""
function movecycle!(st::PersistenceState, (i, j, k))
    empty!(st.σ)
    empty!(st.cycle)
    movepath!(st, i, j)
    movepath!(st, j, k)
    movepath!(st, k, i)
end

"""
    movepath!(st::PersistenceState, u, v)

Move shortest path from `u` to `v` into `st.cycle` and `st.σ`.
"""
function movepath!(st::PersistenceState, u, v)
    while u ≠ v
        u′ = st.parents[v, u]
        push!(st.σ, st.indexmatrix[u, u′])
        push!(st.cycle, u)
        u = u′
    end
end

function isgeodesic(st, Δ, cycle)
    true
end

function processtriangle!(st::PersistenceState, triangle, diameter, showprogress)
    movecycle!(st, triangle)
    # Skip cycles that visit a node more than once.
    length(st.σ) == length(st.cycle) || return
    # Geodesic check
    isgeodesic(st, Δ, cycle) || return

    l = reduce!(st)
    if l ≠ 0
        st.reduced[l] = copy(st.σ)
        if true || length(st.cycle) > 3
            push!(st.generators, copy(st.cycle))
            push!(st.diameters, diameter)
        end
    end
end

function persistence(gc::AbstractGraph, showprogress = false)
    showprogress && println("Calculating intrinsic persistence...")
    st = PersistenceState(gc)

    for (i, (Δ, diam)) in enumerate(st.triangles)
        processtriangle!(st, Δ, diam, showprogress)
    end
    #st
    #filter(x->x≢nothing, map(x->formatresults(gc,x), st.generators))
    showprogress && println("Postprocessing...")
    res = []
    for g in st.generators
        α, m = contract(gc, g)
        if α ≡ nothing
            #println("splat at $(m)!")
        else
            push!(res, (g, α, m))
        end
    end
    res
end

# TODO: optimize
function singlecontract!(α, gc::GeodesicComplex{T}, cycle) where {T}
    changed = false
    ε = zero(T)
    for (i, v) in enumerate(α)
        N = inrange(gc.tree, points(gc, v), gc.radius)
        D = diststocycle(gc, N, cycle)
        ε, j = findmin(vec(maximum(D, dims = 2)))
        changed |= α[i] ≠ N[j]
        α[i] = N[j]
    end
    unique!(α)
    changed, ε
end

function singlecontract_noalloc!(α, gc::GeodesicComplex{T}, cycle) where {T}
    changed = false
    ε = zero(T)
    for (i, u) in enumerate(α)
        N = inrange(gc.tree, points(gc, u), gc.radius)
        # Spaghetti mapreduce avoids allocating distance matrices.
        ε, j = mapreduce(min, enumerate(N)) do (i, v)
            (maximum(evaluate(gc.metric, points(gc, v), landmarks(gc, w))
                     for w in cycle), i)
        end
        changed |= α[i] ≠ N[j]
        α[i] = N[j]
    end
    unique!(α)
    changed, ε
end

function contract(gc::GeodesicComplex{T}, cycle) where {T}
    α = gc.landmarks[cycle]
    ε = zero(T)
    changed = true

    while changed
        changed, ε = singlecontract!(α, gc, cycle)
    end
    length(α) ≥ 3 ? α : nothing, ε
end

diststocycle(gc::GeodesicComplex{T, D}, is, cycle) where {T, D} =
    pairwise(gc.metric,
             reshape(reinterpret(T, points(gc, is)), (D, length(is))),
             reshape(reinterpret(T, landmarks(gc, cycle)), (D, length(cycle))))

#=
struct Cycle
end

struct IntrinsicPersistenceResults{G<:GeodesicComplex}
    complex ::G
    generators ::Vector{...}
    contracted ::Vector{...}
end
=#
