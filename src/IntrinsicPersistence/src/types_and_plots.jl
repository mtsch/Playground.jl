struct IntrinsicGenerator{T}
    triangle ::Triangle{T}
    path     ::Vector{Int}
    lengths  ::Vector{T}
end

IntrinsicGenerator(t::Triangle{T}, path, lens) where T =
    IntrinsicGenerator{T}(t, path, lens)

function Base.show(io::IO, gen::IntrinsicGenerator{T}) where T
    compact = get(io, :compact, false)
    if !compact
        print(io, "IntrinsicGenerator{$T}(")
        showcompact(io, gen.triangle)
        print(io, ", path: $(gen.path), lengths: $(gen.lengths))")
    else
        print(io, "IntrinsicGenerator(")
        showcompact(io, gen.triangle)
        print(io, ", n = $(length(gen.path)), l = $(round(sum(gen.lengths), 5)))")
    end
end

struct IntrinsicPersistenceResult{T, D, I<:IntrinsicMetric{T, D}}
    generators ::Vector{IntrinsicGenerator{T}}
    points     ::Vector{SVector{D, T}}
    subset     ::Vector{Int}
    metric     ::I
end

function Base.show(io::IO, res::IntrinsicPersistenceResult{T}) where T
    println(io, "IntrinsicPersistenceResult{$T}:")
    println(io, "  Number of points: $(length(res.points)) ",
            "($(length(res.subset)) selected)")
    print(io, "  Generators ($(length(res.generators))):")
    for gen in res.generators
        print(io, "\n    ")
        showcompact(io, gen)
    end
end

function Base.convert(::Type{PersistenceBarcode},
                      res::IntrinsicPersistenceResult{T}) where T
    pairs = PersistencePair{T, Void}[]
    for gen in res.generators
        p1 = maximum(gen.lengths)
        p2 = sum(gen.lengths)
        push!(pairs, PersistencePair(p1, p2))
    end
    PersistenceBarcode(PersistencePair{T, Void}[], pairs)
end

# ============================================================================ #
# Plots
# ============================================================================ #
@recipe function f(res::IntrinsicPersistenceResult)
    allpts = res.points
    pts    = allpts[res.subset]
    dim    = length(eltype(pts))
    2 ≤ dim ≤ 3 || throw(ArgumentError("Can't plot $dim-d points!"))

    @series begin
        seriestype := :scatter
        label --> "points"

        xs = getindex.(pts, 1)
        ys = getindex.(pts, 2)
        if dim == 2
            xs, ys
        else
            markersize --> 0.5
            zs = getindex.(pts, 3)
            xs, ys, zs
        end
    end

    @series begin
        seriestype := :scatter
        label --> "all points"
        alpha --> 0.3

        idxs = setdiff(1:size(allpts, 2), res.subset)

        xs = getindex.(allpts[idxs], 1)
        ys = getindex.(allpts[idxs], 2)
        if dim == 2
            xs, ys
        else
            markersize --> 0.5
            zs = getindex.(allpts[idxs], 2)
            xs, ys, zs
        end
    end

    for (i, gen) in enumerate(res.generators)
        @series begin
            l = sum(gen.lengths)

            seriestype := :path
            linewidth --> 5
            label --> "gen$(i); l = $(round(l, 3))"

            path = copy(gen.path)
            push!(path, path[1])

            xs = getindex.(allpts[path], 1)
            ys = getindex.(allpts[path], 2)
            if dim == 2
                xs, ys
            else
                zs = getindex.(allpts[path], 3)
                xs, ys, zs
            end
        end
    end
end
