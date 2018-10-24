using NearestNeighbors
using StaticArrays
using Random
using Plots
using LightGraphs
using SimpleWeightedGraphs
using Distances
using SparseArrays

# some naive samples
function sampler(f)
    function(n)
        res = Vector{SVector{3, Float64}}(undef, n)
        for i in eachindex(res)
            θ, φ = rand(2)
            res[i] = f(θ, φ)
        end
        res
    end
end

torus(R=3, r=1) = sampler((θ, φ) -> @SVector[(R+r*cospi(2θ))cospi(2φ),
                                             (R+r*cospi(2θ))sinpi(2φ),
                                             r*sinpi(2θ)])

sphere(r=1) = sampler((θ, φ) -> @SVector[r*sinpi(θ)cospi(2φ),
                                         r*sinpi(θ)sinpi(2φ),
                                         r*cospi(θ)])

points3d(X; markersize=0.5, kwargs...) =
    scatter3d(getindex.(X, 1), getindex.(X, 2), getindex.(X, 3);
              markersize=markersize, kwargs...)

points2d(X; markersize=1.0, kwargs...) =
    scatter(getindex.(X, 1), getindex.(X, 2);
            markersize=markersize, kwargs...)

points(X; kwargs...) =
    if length(eltype(X)) == 2
        points2d(X; kwargs...)
    else
        points3d(X; kwargs...)
    end

function graphplot(P, G)
    plt = points3d(P)
    E = SVector{3, Float64}[]

    for e in edges(G)
        s = src(e)
        d = dst(e)
        push!(E, P[s])
        push!(E, P[d])
        push!(E, @SVector[NaN,NaN,NaN])
    end

    plot3d!(plt, getindex.(E, 1), getindex.(E, 2), getindex.(E, 3))
end

function uniform_downsample(P, r, kdt = KDTree(X))
    n = length(P)
    kdt = KDTree(P)
    free = trues(n)
    idxs = shuffle(1:n)
    chosen = falses(n)

    for i in idxs
        free[i] || continue
        free[i] = false
        chosen[i] = true
        free[inrange(kdt, P[i], r)] .= false
    end
    (1:n)[chosen]
end

struct GeodesicComplex{T, D, P<:SVector{D, T}, M<:Metric, K<:NNTree{P, M}} <:
                       AbstractSimpleWeightedGraph{Int, T}
    points    ::Vector{P}
    landmarks ::Vector{Int}
    graph     ::SimpleWeightedGraph{Int, T}
    tree      ::K
    cover     ::Vector{Set{Int}}
    radius    ::T
    metric    ::M
end

function GeodesicComplex(pts::AbstractVector{SVector{D, T}}, r;
                         metric=Euclidean(), tree=KDTree) where {D, T}
    kdt = tree(pts, metric)
    landmarks = uniform_downsample(pts, r, kdt)
    cover = [Set(inrange(kdt, pts[landmarks[i]], r)) for i in eachindex(landmarks)]

    I = Int[]
    J = Int[]
    V = T[]
    for i in 1:length(landmarks), j in 1:(i - 1)
        if i ≠ j && !isempty(intersect(cover[i], cover[j]))
            d = evaluate(metric, pts[landmarks[i]], pts[landmarks[j]])
            append!(I, (i, j))
            append!(J, (j, i))
            append!(V, (d, d))
        end
    end
    graph = SimpleWeightedGraph(sparse(I, J, V))

    GeodesicComplex(pts, landmarks, graph, kdt, cover, r, metric)
end

function geodesic_rips(P::Vector{SVector{N, T}}, r) where {N, T}
    S = P[uniform_downsample(P, r)]
    kdt = KDTree(S)
    n = length(S)
    I = Int[]
    J = Int[]
    V = T[]

    for i in 1:n
        for j in inrange(kdt, S[i], 2r)
            v = evaluate(Euclidean(), S[i], S[j])
            append!(I, (i, j))
            append!(J, (j, i))
            append!(V, (v, v))
        end
    end

    S, SimpleWeightedGraph(sparse(I, J, V))
end

function geodesic_witness(P::Vector{SVector{N, T}}, r) where {N, T}
    kdt = KDTree(P)
    chosen = uniform_downsample(P, r)
    n = length(chosen)
    I = Int[]
    J = Int[]
    V = T[]

    for i in chosen
        for j in chosen
            if !isempty(intersect(inrange(kdt, P[i], r), inrange(kdt, P[j], r)))
                v = evaluate(Euclidean(), P[i], P[j])
                append!(I, (i, j))
                append!(J, (j, i))
                append!(V, (v, v))
            end
        end
    end

    P[chosen], SimpleWeightedGraph(sparse(I, J, V))
end


# generate some plots
function generate(n; r = 0.2, out="img", ext=".pdf")
    S2 = sphere()(n)
    T2 = torus()(n)
    G  = [@SVector[randn(), randn()] for _ in 1:n]

    for (X, name) in zip([S2, T2, G], ["S^2", "T^2", "gauss"])
        xs = getindex.(X, 1)
        ys = getindex.(X, 2)
        if length(eltype(X)) > 2
            zs = getindex.(X, 3)
        else
            zs = 0
        end
        xlim = extrema(xs)
        ylim = extrema(xs)
        zlim = extrema(zs)

        plt = points(X, #title="$name n=$(length(X))",
                     xlim=xlim, ylim=ylim, zlim=zlim,
                     aspect_ratio = 1,
                     legend = false)
        savefig(plt, joinpath(out, "$(name)_full$ext"))

        Y = uniform_downsample(X, r)
        plt = points(Y, #title="$name uniform r=$r n=$(length(Y))",
                     xlim=xlim, ylim=ylim, zlim=zlim,
                     aspect_ratio = 1,
                     legend = false)
        savefig(plt, joinpath(out, "$(name)_uniform$ext"))

        Y = weighted_downsample(X, r)
        plt = points(Y, #title="$name weighted r=$r n=$(length(Y))",
                     xlim=xlim, ylim=ylim, zlim=zlim,
                     aspect_ratio = 1,
                     legend = false)
        savefig(plt, joinpath(out, "$(name)_weighted$ext"))
    end
end

normalize.(@SVector rand(3) for _ in 1:1000) |> points


function datasets()

    T3 = torus()(1_000_000)
    g3 = GeodesicComplex(T3, 1.)

    T4 = vcat(T3[1:500_000], map(x -> x.+2.5, T3[500_0001:end]))
    g4w = GeodesicComplex(T4, .5)
    g4r = GeodesicComplex(T4, .5, witness = false)
end
