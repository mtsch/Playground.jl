# ============================================================================ #
# DictSet
# ============================================================================ #
# TODO: use IntSet? - benchmark!
"""
Basic dict-like thing that holds multiple values for each key in a set.
"""
struct DictSet{K}
    data::Dict{K, IntSet}
end

DictSet{K}() where K = DictSet{K}(Dict{K, IntSet}())

function Base.getindex(ds::DictSet{K}, k::K) where K
    if haskey(ds.data, k)
        ds.data[k]
    else
        IntSet()
    end
end

function Base.insert!(ds::DictSet{K}, k::K, v::Int) where K
    if haskey(ds.data, k)
        push!(ds.data[k], v)
    else
        ds.data[k] = IntSet(v)
    end
    v
end

function Base.delete!(ds::DictSet{K}, k::K) where K
    delete!(ds.data, k)
    ds
end

function Base.delete!(ds::DictSet{K}, k::K, v::Int) where K
    if haskey(ds.data, k)
        s = ds.data[k]
        if length(s) > 1
            delete!(s, v)
        else
            delete!(ds.data, k)
        end
    end
    ds
end

function Base.delete!(ds::DictSet{K}, vals::IntSet) where K
    for (k, v) in ds.data
        setdiff!(v, vals)
        if isempty(v)
            delete!(ds.data, k)
        end
    end
end

Base.keys(ds::DictSet) = keys(ds.data)
Base.values(ds::DictSet) = union(values(ds.data)...)
Base.haskey(ds::DictSet) = haskey(ds.data)
Base.isempty(ds::DictSet) = isempty(ds.data)

# ============================================================================ #
# THE SAL TYPE
# ============================================================================ #
mutable struct SAL{T}
    nodelist::Vector{DictSet{Int}}

    # TODO: store data
    # reverse dict for maximal simplices?
    lbl_vals::IntSet
    lbl_holes::IntSet
end

SAL{T}() where {T} = SAL{T}(DictSet{Int}[], IntSet(), IntSet())
SAL() = SAL{Float64}()

Base.show(io::IO, sal::SAL) = print(io, "SAL with $(length(sal)) nodes")

# Pretty printing.
function Base.print(io::IO, sal::SAL)
    println("SAL:")
    for i in eachindex(sal)
        print("$(lpad(i, 5)): [")
        for k in sort(collect(keys(sal.nodelist[i].data)))
            for v in sort(collect(sal.nodelist[i].data[k]))
                print("($i, $(iszero(k) ? "φ" : k), $v)")
            end
        end
        println("]")
    end
end

# ============================================================================ #
# GETTERS, SETTERS, UTILS
# ============================================================================ #
function getlabel!(sal::SAL)
    if isempty(sal.lbl_holes)
        lbl = length(sal.lbl_vals) + 1
    else
        lbl = first(sal.lbl_holes)
        delete!(sal.lbl_holes, lbl)
    end
    push!(sal.lbl_vals, lbl)
    lbl
end

Base.length(sal::SAL) = length(sal.nodelist)
Base.eachindex(sal::SAL) = 1:length(sal.nodelist)
Base.isempty(sal::SAL) = all(isempty.(sal.nodelist))

function insnode!(sal::SAL, cur::Int, nxt::Int, lbl::Int)
    while length(sal) < cur
        push!(sal.nodelist, DictSet{Int}())
    end

    insert!(sal.nodelist[cur], nxt, lbl)
end

function labels(sal::SAL, from::Int, to::Int)
    sal.nodelist[from][to]
end

function alllabels(sal::SAL)
    sal.lbl_values
end

function maximalsimplices(sal::SAL)
    maxsxs = Dict{Int, Vector{Int}}()
    for l in sal.lbl_vals
        maxsxs[l] = []
    end
    for i in eachindex(sal)
        for v in values(sal.nodelist[i])
            push!(maxsxs[v], i)
        end
    end
    maxsxs
end

"""
Check if σ1 is a face of σ2. Assumes the simplices are properly formatted.
"""
function isface(σ1, σ2)
    i = 1; j = 1
    while i <= length(σ1) && j <= length(σ2)
        if σ1[i] > σ2[j]
            j += 1
        elseif σ1[i] == σ2[j]
            i += 1; j += 1
        else
            return false
        end
    end
    i == length(σ1) + 1
end

# ============================================================================ #
# OPERATIONS
# ============================================================================ #
# A σ is in the complex if it is contained in one of the maximal sx's.
function Base.in(σ, sal::SAL)
    σ = check_simplex(σ)
    !isempty(inmaximals(sal, σ))
end

"""
List maximal simplices that contain σ.
"""
function inmaximals(sal::SAL, σ)
    last(σ) > length(sal) && return IntSet()

    s = copy(labels(sal, σ[1], σ[2]))
    for l in zip(σ[2:end-1], σ[3:end])
        intersect!(s, labels(sal, l...))
        isempty(s) && return s
    end
    s
end

function Base.push!(sal::SAL, σ)
    # TODO: remove faces of simplex!
    # Find maximal simplices that are contained in σ and remove them
    σ = check_simplex(σ)
    σ in sal && return σ

    lbl = getlabel!(sal)

    for (i, li) in enumerate(σ), (j, lj) in enumerate(σ[i+1:end])
        insnode!(sal, li, lj, lbl)
    end
    insnode!(sal, σ[end], 0, lbl)
    σ
end

function deletemaximals!(sal::SAL, lbls::IntSet)
    for i in eachindex(sal)
        delete!(sal.nodelist[i], lbls)
    end
    union!(sal.lbl_holes, lbls)
    setdiff!(sal.lbl_vals, lbls)
end

"""

"""
function facetswithout(σ, τ)
end

function Base.delete!(sal::SAL, σ)
    # TODO: multiple at once!
    # insert facets that do not contain σ
    deleted_lbls = IntSet()
    deleted_sxs  = Vector{Int}[]
    for (lbl, τ) in maximalsimplices(sal)
        if isface(τ, σ)
            push!(delted_lbls, lbl)
            push!(delted_sxs, τ)
        end
    end
    deletemaximals!(sal, deleted_lbls)
    for τ in delted_sxs
        for υ in facetswithout(τ, σ)
            insert!(sal, υ)
        end
    end

    # vsi k vsebujejo σ = greš po drevesu kot σ in delaš intersectione
    for (i, li) in enumerate(σ)
    end
end
