module GeodesicComplexes
export GeodesicComplex, npoints, nlandmarks, points, landmarks, nonlandmarks

using Distances
using LightGraphs
using NearestNeighbors
using RecipesBase
using SimpleWeightedGraphs
using StaticArrays

using Random
using SparseArrays

struct GeodesicComplex{T, D, P<:SVector{D, T}, M<:Metric, K<:NNTree{P, M}} <:
                       AbstractSimpleWeightedGraph{Int, T}
    points    ::Vector{P}
    landmarks ::Vector{Int}
    graph     ::SimpleWeightedGraph{Int, T}
    tree      ::K # is this needed?
    cover     ::Vector{Set{Int}}
    radius    ::T
    metric    ::M
    witness   ::Bool
end

# TODO:
# Mogoče je treba r ---> 2r, ker drugače lahko dobiš trikotnike, ki ne morejo imet
# presečišča, ne glede na gostoto točk. Je 2r vedno dovolj?
function GeodesicComplex(pts::AbstractVector{SVector{D, T}}, r;
                         metric = Euclidean(), tree = KDTree,
                         witness = true) where {D, T}
    kdt = tree(pts, metric)
    landmarks, cover = getcover(pts, r, kdt)
    n = length(landmarks)

    I = Int[]
    J = Int[]
    V = T[]
    for i in 1:n, j in 1:(i - 1)
        d = evaluate(metric, pts[landmarks[i]], pts[landmarks[j]])
        if (witness && !isempty(intersect(cover[i], cover[j]))) ||
           (!witness && d < 2r)
            append!(I, (i, j))
            append!(J, (j, i))
            append!(V, (d, d))
        end
    end
    graph = SimpleWeightedGraph(sparse(I, J, V, n, n))

    GeodesicComplex(pts, landmarks, graph, kdt, cover, T(r), metric, witness)
end

function getcover(P, r, kdt = KDTree(P))
    covered = falses(length(P))
    idxs = shuffle(eachindex(P))
    landmarks = Int[]
    cover = Set{Int}[]

    for i in idxs
        covered[i] && continue
        covered[i] = true
        # names!
        inball = inrange(kdt, P[i], r)
        covered[inball] .= true
        push!(landmarks, i)
        push!(cover, Set(inball))
    end
    landmarks, cover
end

# LightGraphs interface
for f in [:edges, :eltype, :ne, :nv, :vertices, :is_directed, :weights]
    @eval LightGraphs.$f(gc::GeodesicComplex) = $f(gc.graph)
end
LightGraphs.is_directed(::Type{GeodesicComplex}) = false
LightGraphs.is_directed(::Type{GeodesicComplex{T,D,P,M,K}}) where {T,D,P,M,K} = false
LightGraphs.has_edge(gc::GeodesicComplex, u::Integer, v::Integer) = has_edge(gc.graph, u, v)
LightGraphs.has_vertex(gc::GeodesicComplex, v::Integer) = has_edge(gc.graph, v)
LightGraphs.inneighbors(gc::GeodesicComplex, v::Integer) = inneighbors(gc.graph, v)
LightGraphs.outneighbors(gc::GeodesicComplex, v::Integer) = outneighbors(gc.graph, v)

function Base.show(io::IO, gc::GeodesicComplex{T, D}) where {T, D}
    print(io, "GeodesicComplex{$T, $D} with $(npoints(gc)) points, " *
          "$(nlandmarks(gc)) landmarks and $(ne(gc)) edges")
end

npoints(gc::GeodesicComplex) =
    length(gc.points)
nlandmarks(gc::GeodesicComplex) =
    length(gc.landmarks)
points(gc::GeodesicComplex, idxs=1:npoints(gc)) =
    gc.points[idxs]
landmarks(gc::GeodesicComplex, idxs=1:nlandmarks(gc)) =
    gc.points[gc.landmarks[idxs]]
nonlandmarks(gc::GeodesicComplex, idxs=1:npoints(gc)-nlandmarks(gc)) =
    gc.points[setdiff(1:npoints(gc), gc.landmarks)[idxs]]

# Plots recipe
function getxyz(pts)
    xs = get.(pts, 1, 0.0)
    ys = get.(pts, 2, 0.0)
    zs = get.(pts, 3, 0.0)
    if all(iszero, zs)
        xs, ys
    else
        xs, ys, zs
    end
end

@recipe function plot(gc::GeodesicComplex{T, D}, cycles = [];
                      only_landmarks = true, graph = true) where {T, D}
    # edges
    if graph
        @series begin
            label := "edges"

            edgepoints = SVector{D, Float64}[]
            for e in edges(gc)
                s = landmarks(gc, src(e))
                d = landmarks(gc, dst(e))
                append!(edgepoints, (s, d, @SVector fill(NaN, D)))
            end

            getxyz(edgepoints)
        end
    end

    # landmarks
    @series begin
        markersize --> 1.0
        seriestype --> :scatter
        label := "landmarks"

        getxyz(landmarks(gc))
    end

    # others
    if !only_landmarks
        @series begin
            markersize --> 0.5
            seriestype --> :scatter
            label := "others"

            getxyz(nonlandmarks(gc))
        end
    end

    # cycles
    for (i, c) in enumerate(cycles)
        @series begin
            label := "cycle $i"
            linewidth := 5
            edgepoints = vcat(landmarks(gc, c), [landmarks(gc, c[1])])
            getxyz(edgepoints)
        end
    end
end

end
