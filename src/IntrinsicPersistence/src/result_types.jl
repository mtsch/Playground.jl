# ============================================================================ #
# RESULT TYPES
# ============================================================================ #
struct IntrinsicGenerator{T}
    triangle ::Triangle{T}
    path     ::Vector{Int}
    lengths  ::Vector{T}
end

#IntrinsicGenerator(t::Triangle{T}, path, lens) where T =
#    IntrinsicGenerator{T}(t, path, lens)

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

Base.length(gen::IntrinsicGenerator) = sum(gen.lengths)

struct IntrinsicPersistenceResult{T, D, G<:AbstractGraph}
    generators ::Vector{IntrinsicGenerator{T}}
    graph      ::G
    data       ::D

function IntrinsicPersistenceResult(gens  ::Vector{IntrinsicGenerator{T}},
                                    graph ::G,
                                    data  ::D = nothing
                                    ) where {T, G<:AbstractGraph, D}

    new{T, D, G}(gens, graph, data)
end
end


function Base.show(io::IO, res::IntrinsicPersistenceResult{T, D}) where {T, D}
    println(io, "IntrinsicPersistenceResult{$T,$D}:")
    println(io, "  Graph: $(res.graph)")
    print(io,   "  Generators ($(length(res.generators))):")
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

# Associate some data with the result. This is used to add points to a result
# that was calculated from a spanner graph.
add_data(res::IntrinsicPersistenceResult{<:Any, Void, <:Any}, data) =
    IntrinsicPersistenceResult(res.generators, res.graph, data)

# ============================================================================ #
# Plots
# ============================================================================ #
# Get x, y, and z coordinates from matrix or arrav of SVectors.
function get_xyz(pts::AbstractMatrix{<:Any}, i=Colon())
    dim = size(pts, 1)
    dim > 3 && warn("Data is in $(dim)d. Only using the first three coordinates.",
                    once = true)
    if dim == 1
        x,y,z = pts[1, i], zeros(size(pts, 2)), zeros(size(pts, 2))
    elseif dim == 2
        x,y,z = pts[1, i], pts[2, i], zeros(size(pts, 2))
    else
        x,y,z = pts[1, i], pts[2, i], pts[3, i]
    end
    Float64.(x), Float64.(y), Float64.(z)
end

function get_xyz(pts::AbstractVector{<:StaticArray}, i=Colon())
    dim = length(eltype(pts))
    dim > 3 && warn("Data is in $(dim)d. Only using the first three coordinates.",
                    once = true)
    pts = pts[i]
    if dim == 1
        x,y,z = getindex.(pts, 1), zeros(length(pts)), zeros(length(pts))
    elseif dim == 2
        x,y,z = getindex.(pts, 1), getindex.(pts, 2), zeros(length(pts))
    else
        x,y,z = getindex.(pts, 1), getindex.(pts, 2), getindex.(pts, 3)
    end
    Float64.(x), Float64.(y), Float64.(z)
end

function get_xyz(::T) where T
    throw(ArgumentError("Can't plot $T!"))
end

@recipe function f(res::IntrinsicPersistenceResult{<:Any}, graph = 0)
    aspect_ratio := 1

    graph ≥ 2 && @series begin
        xs = Float64[]; ys = Float64[]; zs = Float64[]
        for Δ in equilaterals(res.graph.weights, tol = Inf)
            i,j,k = Δ.vertices
            xyz = get_xyz(res.data, [i,j,k,i])
            append!(xs, xyz[1]); push!(xs, NaN)
            append!(ys, xyz[2]); push!(ys, NaN)
            append!(zs, xyz[3]); push!(zs, NaN)
        end

        seriestype := :shape
        alpha := 0.1
        label := "2-skeleton"

        if all(z -> z == 0 || z == NaN, zs)
            xs, ys
        else
            xs, ys, zs
        end
    end

    graph ≥ 1 && @series begin
        xs = Float64[]; ys = Float64[]; zs = Float64[]
        for e in edges(res.graph)
            p1 = src(e)
            p2 = dst(e)
            xyz = get_xyz(res.data, [p1, p2])
            append!(xs, xyz[1]); push!(xs, NaN)
            append!(ys, xyz[2]); push!(ys, NaN)
            append!(zs, xyz[3]); push!(zs, NaN)
        end

        markersize := :none
        seriestype := :path
        linewidth --> 0.5
        if graph ≥ 2
            color := :black
        end
        label := "1-skeleton"

        if all(z -> z == 0 || z == NaN, zs)
            xs, ys
        else
            xs, ys, zs
        end
    end

    @series begin
        seriestype := :scatter
        label --> :points
        markersize --> 0.5

        xs, ys, zs = get_xyz(res.data)
        if all(z -> z == 0 || z == NaN, zs)
            xs, ys
        else
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

            xs, ys, zs = get_xyz(res.data, path)
            if all(z -> z == 0 || z == NaN, zs)
                xs, ys
            else
                xs, ys, zs
            end
        end
    end
end
