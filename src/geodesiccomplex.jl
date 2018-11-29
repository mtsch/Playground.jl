module GeodesicComplexes
export GeodesicComplex, n_points, n_landmarks, points, landmark_points, nonlandmark_points

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

function GeodesicComplex(pts::AbstractVector{SVector{D, T}}, r;
                         metric=Euclidean(), tree=KDTree, witness=true) where {D, T}
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
    print(io, "GeodesicComplex{$T, $D} with $(n_points(gc)) points, " *
          "$(n_landmarks(gc)) landmarks and $(ne(gc)) edges")
end

n_points(gc) = length(points(gc))
n_landmarks(gc) = length(gc.landmarks)
points(gc) = gc.points
landmark_points(gc::GeodesicComplex) = gc.points[gc.landmarks]
nonlandmark_points(gc::GeodesicComplex) = gc.points[setdiff(1:n_points(gc), gc.landmarks)]

# Plots recipe
function getxyz(pts::AbstractVector{SVector{D, T}}) where {D, T}
    xs = getindex.(pts, 1)
    ys = getindex.(pts, 2)
    if D ≥ 3
        zs = getindex.(pts, 3)
        xs, ys, zs
    else
        xs, ys
    end
end

@recipe function plot(gc::GeodesicComplex{T, D}, cycles = []; only_landmarks = true) where {T, D}
    landmarks = landmark_points(gc)

    # edges
    @series begin
        label := "edges"

        edgepoints = SVector{D, Float64}[]
        for e in edges(gc)
            s = landmarks[src(e)]
            d = landmarks[dst(e)]
            append!(edgepoints, (s, d, @SVector fill(NaN, D)))
        end

        getxyz(edgepoints)
    end

    # landmarks
    @series begin
        markersize --> 1.0
        seriestype --> :scatter
        label := "landmarks"

        getxyz(landmark_points(gc))
    end

    # others
    if !only_landmarks
        @series begin
            markersize --> 0.5
            seriestype --> :scatter
            label := "others"
            alpha --> 0.5

            getxyz(nonlandmark_points(gc))
        end
    end

    # cycles
    for (i, c) in enumerate(cycles)
        @series begin
            label := "cycle $i"
            linewidth := 5
            edgepoints = vcat(landmarks[c], [landmarks[c[1]]])
            getxyz(edgepoints)
        end
    end
end

end
