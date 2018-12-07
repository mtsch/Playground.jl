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

# TODO: islocallyminimal
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
    error()
end

function processtriangle!(st::PersistenceState, triangle, diameter, showprogress)
    movecycle!(st, triangle)
    # Skip cycles that visit a node more than once.
    length(st.σ) == length(st.cycle) || return
    # Geodesic check

    # showprogress && print(triangle, " ↦ ", st.cycle, ", ", collect(st.σ))
    l = reduce!(st)
    if l ≠ 0
        st.reduced[l] = copy(st.σ)
        if true || length(st.cycle) > 3
            push!(st.generators, copy(st.cycle))
            push!(st.diameters, diameter)
        end
        #showprogress && println(" survived as ", collect(st.σ), ".")
    else
        #showprogress && println(" is kill.")
    end
end

function persistence(gc::AbstractGraph, showprogress = false)
    showprogress && println("Calculating intrinsic persistence...")
    st = PersistenceState(gc)

    @time for (i, (Δ, diam)) in enumerate(st.triangles)
        processtriangle!(st, Δ, diam, showprogress)
    end
    #st
    #filter(x->x≢nothing, map(x->formatresults(gc,x), st.generators))
    showprogress && println("Postprocessing...")
    res = []
    @time for g in st.generators
        α, nit, m = contract(gc, g)
        if α ≡ nothing
            #println("splat at $(m)!")
            if m>2gc.radius
                println("ujjjjjjjj, $m")
            end
        else
            push!(res, (g, centroid(gc, g)[1], α, nit, m))
        end
    end
    res
    #st
end

function contract(gc, cycle)
    β = gc.landmarks[cycle]
    α = Int[]
    m = 0.0

    nit = 0
    while α ≠ β
        resize!(α, length(β))
        copyto!(α, β)
        for (i, v) in enumerate(α)
            N = inrange(gc.tree, points(gc, v), gc.radius)
            D = diststocycle(gc, N, cycle)
            m, j = findmin(vec(maximum(D, dims = 2)))
            β[i] = N[j]
        end
        unique!(β)
        nit += 1
        if length(β) ≤ 2
            return nothing, nit, m
        end
    end
    β, nit, m
end

_pairwise(gc::GeodesicComplex{T, D}, is, js) where {T, D} =
    pairwise(gc.metric,
             reshape(reinterpret(T, points(gc, is)), (D, length(is))),
             reshape(reinterpret(T, points(gc, js)), (D, length(js))))

diststocycle(gc::GeodesicComplex{T, D}, is, cycle) where {T, D} =
    pairwise(gc.metric,
             reshape(reinterpret(T, points(gc, is)), (D, length(is))),
             reshape(reinterpret(T, landmarks(gc, cycle)), (D, length(cycle))))

function centroid(gc, cycle)
    r = gc.radius
    tree = gc.tree
    metric = gc.metric

    N = collect(mapreduce(i -> gc.cover[i], union, cycle))
    D = diststocycle(gc, N, cycle)
    d, i = findmin(vec(maximum(D, dims = 2)))
    p = N[i]
    d′ = typemax(d)
    while d′ > d
        d′ = d
        N = inrange(tree, points(gc, p), r)
        D = diststocycle(gc, N, cycle)
        d, i = findmin(vec(maximum(D, dims = 2)))
        p = N[i]
    end
    d, p
end

tomatrix(arr::AbstractArray{SVector{N, T}}) where {N, T} =
    reshape(reinterpret(T, arr), (N, length(arr)))

#=
formatresults(g::AbstractGraph, cycle) = cycle

function formatresults(gc::GeodesicComplex, cycle)
    # Tole je narobe, ampak ideja je:
    # najdeš točko, ki je najbližje točkam iz cikla
    # najdaljša razdalja do nje je radij smrti
    # POZOR: točka smrti ni nujno v okolici cikla!
    surrounding = mapreduce(i -> gc.cover[i], union, cycle)
    all_pts = tomatrix(points(gc)[collect(surrounding)])
    lnd_pts = tomatrix(landmark_points(gc)[cycle])
    dists = pairwise(gc.metric, all_pts, lnd_pts)
    maxs = findmin(maximum(dists, dims = 1))
    death, ideath = findmin()

    if death < gc.radius #*2#?
        nothing
    else
        cycle, death
    end
end
=#

function cleanup(gc::GeodesicComplex, st::PersistenceState{T};
                 removetriangles = false) where {T}
    generators = Vector{Int}[]
    diameters = T[]
    for g in st.generators
        if !removetriangles || length(g) > 3
            if isempty(mapreduce(i -> gc.cover[i], intersect, g))
                push!(generators, g)
            end
        end
    end
    generators
end
