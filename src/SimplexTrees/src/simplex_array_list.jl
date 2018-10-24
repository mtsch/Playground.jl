# ============================================================================ #
# CORE TYPES
# ============================================================================ #
"""
Basic dict-like thing that holds multiple values for each key in a set.
"""
struct DictSet{K, V}
    data::Dict{K, Set{V}}
end

DictSet{K, V}() where {K, V} = DictSet{K, V}(Dict{K, V}())

function Base.getindex(ds::DictSet{K, V}, k::K) where {K, V}
    if haskey(ds.data, k)
        ds.data[k]
    else
        Set{V}()
    end
end

function Base.insert!(ds::DictSet{K, V}, k::K, v::V) where {K, V}
    if haskey(ds.data, k)
        push!(ds.data[k], v)
    else
        ds.data[k] = Set{V}(v)
    end
    v
end

function Base.delete!(ds::DictSet{K, V}, k::K) where {K, V}
    delete!(ds.data, k)
    ds
end

function Base.delete!(ds::DictSet{K, V}, k::K, v::V) where {K, V}
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

"""

"""
struct SALNode{L, T}
    next_lbl::L
    mxsx_lbl::L

    data::Nullable{T}
end

function Base.isless(n1::SALNode, n2::SALNode)
    for f in [:next_lbl, :mxsx_lbl]
        getfield(n1, f) < getfield(n2, f) && return true
        getfield(n1, f) > getfield(n2, f) && return false
    end
    n1.data < n2.data
end

SALNode{L, T}(next_lbl::L, mxsx_lbl::L) where {L, T} =
    SALNode(next_lbl, mxsx_lbl, Nullable{T}())

Base.show(io::IO, nd::SALNode) =
    print(io, "(" * (iszero(nd.next_lbl) ? "φ" : "$(nd.next_lbl)") *
          ", $(nd.mxsx_lbl)" * (isnull(nd.data) ? "" : ", $(nd.data)") * ")")

"""

"""
mutable struct SAL{L, T}
    nodes::Vector{SortedSet{SALNode{L, T}}}

    last_label::L
end

function SAL{L, T}() where {L, T}
    SAL(SortedSet{SALNode{L, T}}[], 0)
end
SAL(::Type{L}, ::Type{T}) where {L, T} = SAL{L, T}()
SAL() = SAL(Int, Float64)

Base.show(io::IO, sal::SAL) =
    print("SAL with $(n_vertices(sal)) vertices.")

# Pretty print, showing structure.
function Base.print(io::IO, sal::SAL)
    println("SAL:")
    for i in 1:n_vertices(sal)
        print("$(lpad(i, 5)): [")
        nds = collect(sal.nodes[i])
        for el in nds[1:end-1]
            print("$el, ")
        end
        println("$(nds[end])]")
    end
end

# ============================================================================ #
# GETTERS, SETTERS, UTILS
# ============================================================================ #
get_label!(sal::SAL) = sal.last_label += 1

n_vertices(sal::SAL) = length(sal.nodes)

"""
Create a new node and insert it into the SAL.
"""
function insert_node!(sal::SAL{L, T}, ndlbl, nxlbl, sxlbl) where {L, T}
    while n_vertices(sal) < ndlbl
        push!(sal.nodes, SortedSet{SALNode{L, T}}())
    end

    push!(sal.nodes[ndlbl], SALNode{L, T}(nxlbl, sxlbl))
end

# ============================================================================ #
# OPERATIONS
# ============================================================================ #
function insert_simplex!(sal::SAL, simplex)
    # TODO: rename to Base.push?
    # TODO: remove faces of simplex!
    # Find maximal simplices that are contained in σ and remove them

    check_simplex(simplex)
    n = length(simplex)

    sxlbl = get_label!(sal)

    for (i, li) in enumerate(simplex), (j, lj) in enumerate(simplex[i+1:end])
        insert_node!(sal, li, lj, sxlbl)
    end
    insert_node!(sal, simplex[end], 0, sxlbl)
end

function Base.in(simplex, sal::SAL)
    simplex[end] > n_vertices(sal) && return false
    check_simplex(simplex)

    # TODO: filter is inefficient!
    findmatches(node, next) =
        IntSet(map(x -> x.mxsx_lbl,
                   Iterators.filter(x -> x.next_lbl == next, sal.nodes[node])))

    s = findmatches(simplex[1], simplex[2])
    for l in zip(simplex[2:end-1], simplex[3:end])
        intersect!(s, findmatches(l...))
        isempty(s) && return false
    end
    true
end

function remove_simplex!(sal::SAL, simplex)
    # find maximal simplices which contain σ
    # remove them
    # insert facets that do not contain σ
end
