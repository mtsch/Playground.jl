# ============================================================================ #
# The SimplicialComplex type
# ============================================================================ #
"""
  SimplicialComplex{T}

The simplicial complex. Holds data of type `T` (not implemented.)

Supported operations:

* push!(sxcx, σ)   [√]
* σ in sxcx        [√]
* delete!(sxcx, σ) [Χ]

"""
mutable struct SimplicialComplex{T}
    nodelist::Vector{DictSet{Int, Int}}

    # TODO: store data
    # reverse dict for maximal simplices?
    lbl_values::IntSet
    lbl_holes::IntSet
end

SimplicialComplex{T}() where T =
    SimplicialComplex{T}(DictSet{Int, Int}[], IntSet(), IntSet())
SimplicialComplex() =
    SimplicialComplex{Float64}()

# ---------------------------------------------------------------------------- #
# Base overloads.

Base.show(io::IO, sxcx::SimplicialComplex) =
    print(io, "SimplicialComplex with $(length(sxcx)) nodes")

# Pretty printing.
function Base.print(io::IO, sxcx::SimplicialComplex)
    println(io, "SimplicialComplex:")
    for i in eachindex(sxcx)
        print(io, "$(lpad(i, 5)): [")
        for k in sort(collect(keys(sxcx.nodelist[i].data)))
            for v in sort(collect(sxcx.nodelist[i].data[k]))
                print(io, "($i, $(iszero(k) ? "φ" : k), $v)")
            end
        end
        println(io, "]")
    end
end

Base.length(sxcx::SimplicialComplex) = length(sxcx.nodelist)
Base.eachindex(sxcx::SimplicialComplex) = 1:length(sxcx.nodelist)
Base.isempty(sxcx::SimplicialComplex) = all(isempty.(sxcx.nodelist))

# ============================================================================ #
# Operation helpers.
# ============================================================================ #
"""
  newlabel!(sxcx)

Create a new maximal simplex label and return it.
"""
function newlabel!(sxcx::SimplicialComplex)
    if isempty(sxcx.lbl_holes)
        lbl = length(sxcx.lbl_values) + 1
    else
        lbl = first(sxcx.lbl_holes)
        delete!(sxcx.lbl_holes, lbl)
    end
    push!(sxcx.lbl_values, lbl)
    lbl
end

"""
  insnode!(sxcx, cur, nxt, lbl)

Insert a node that represents a link between `cur` and `nxt` in the maximal
simplex with label `lbl`.
"""
function insnode!(sxcx::SimplicialComplex, cur::Int, nxt::Int, lbl::Int)
    while length(sxcx) < cur
        push!(sxcx.nodelist, DictSet{Int, Int}())
    end

    insert!(sxcx.nodelist[cur], nxt, lbl)
end

"""
  labels(sxcx, from, to)

Return the labels of maximal simplices that contain the link between `from` and
`to`.
"""
function labels(sxcx::SimplicialComplex, from::Int, to::Int)
    sxcx.nodelist[from][to]
end

"""
  alllabels(sxcx)

Return all labels of maximal simplices in `sxcx`.
"""
function alllabels(sxcx::SimplicialComplex)
    sxcx.lbl_values
end

"""
  maximalsimplices(sxcx)

Return all maximal simplices in `sxcx` in a Dict{Int, Vector{Int}}.

TODO: This is slow, make an iterator or avoid doing this entirely.
"""
function maximalsimplices(sxcx::SimplicialComplex)
    maxsxs = Dict{Int, Vector{Int}}()
    for l in sxcx.lbl_values
        maxsxs[l] = []
    end
    for i in eachindex(sxcx)
        for v in values(sxcx.nodelist[i])
            push!(maxsxs[v], i)
        end
    end
    maxsxs
end

# ============================================================================ #
# Operations
# ============================================================================ #
# A σ is in the complex if it is contained in one of the maximal sx's.
function Base.in(σ, sxcx::SimplicialComplex)
    σ = validatesx(σ)
    !isempty(maximalcofaces(sxcx, σ))
end

"""
  maximalcofaces(sxcx, σ)

List the labels of maximal simplices that are cofaces of `σ`.
"""
function maximalcofaces(sxcx::SimplicialComplex, σ)
    last(σ) > length(sxcx) && return Set{Int}()

    s = labels(sxcx, σ[1], σ[2])
    for l in zip(σ[2:end-1], σ[3:end])
        s = intersect(s, labels(sxcx, l...))
        isempty(s) && return s
    end
    s
end

"""
  maximalfaces(sxcx, σ)

List the labels of maximal simplices that are faces of `σ`.
"""
function maximalfaces(sxcx::SimplicialComplex, σ)
    # Find labels of maximal simplices that have a node that is not in σ and
    # setdiff them from all labels.
    s = Set{Int}()
    for i in setdiff(eachindex(sxcx), σ)
        union!(s, values(sxcx.nodelist[i]))
    end

    maxsxs = Set(sxcx.lbl_values)
    setdiff(maxsxs, s)
end

"""
  deletemaximals!(sxcx, lbls)

Delete maximal simplices with labels in lbls.
"""
function deletemaximals!(sxcx::SimplicialComplex, lbls::Set{Int})
    for i in eachindex(sxcx)
        delete!(sxcx.nodelist[i], lbls)
    end
    union!(sxcx.lbl_holes, lbls)
    setdiff!(sxcx.lbl_values, lbls)
end

# DONE.
function Base.push!(sxcx::SimplicialComplex, σ)
    σ = validatesx(σ)
    σ in sxcx && return σ

    # Remove maximal simplices that are faces of σ.
    deletemaximals!(sxcx, maximalfaces(sxcx, σ))

    lbl = newlabel!(sxcx)
    for (i, li) in enumerate(σ), (j, lj) in enumerate(σ[i+1:end])
        insnode!(sxcx, li, lj, lbl)
    end
    insnode!(sxcx, σ[end], 0, lbl)
    σ
end

"""

"""
function facetswithout(σ, τ)
end

function Base.delete!(sxcx::SimplicialComplex, σ)
    maxlbls = maximalcofaces(sxcx, σ)
    # make dim(σ) copies of each of the maximal simplices and in the i-th copy
    # remove nodes with label (σ[i], x, ℓ)

    # Torej v bistvu iz vsakega daš eno točko
end
