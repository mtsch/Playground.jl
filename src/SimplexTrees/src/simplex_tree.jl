# ============================================================================ #
# CORE TYPES
# ============================================================================ #
"""
    AbstractNode{L<:Integer, T}

Can be a `SXTree` or a `SXTNode`. `L` is the type of label and `T` is the type
of data.
"""
abstract type AbstractNode{L<:Integer, T} end

"""
    SXTNode{L, T} <: AbstractNode{L, T}

A node in a simplex tree. The `data::T` field is nullable and empty by default.

Fields in constructor:

* `label::L`
* `depth::L`
* `parent::AbstractNode{L, T}`
* `data::Nullable{T}`
"""
struct SXTNode{L, T} <: AbstractNode{L, T}
    label::L
    depth::L

    data::Nullable{T}

    parent::AbstractNode{L, T}
    children::Dict{L, SXTNode{L, T}}

    next::Base.RefValue{SXTNode{L, T}}
    prev::Base.RefValue{SXTNode{L, T}}
end

function SXTNode(l::L, d::L,
              p::AbstractNode{L, T}, data=Nullable{T}()) where {L<:Integer, T}

    v = SXTNode(l, d, data, p, Dict{L, SXTNode{L, T}}(),
             Base.RefValue{SXTNode{L, T}}(), Base.RefValue{SXTNode{L, T}}())
    v.next[] = v
    v
end

Base.show(io::IO, nd::SXTNode) =
    print(io, "SXTNode{label=$(label(nd)), depth=$(depth(nd))",
          !isnull(nd.data) ? ", data=$(get(nd.data))" : "", "}")


"""
SXTree{L, T} <: AbstractNode{L, T}

The simplex tree.

TODO
"""
struct SXTree{L, T} <: AbstractNode{L, T}
    children::Vector{SXTNode{L, T}}
    cousin_list::Vector{Dict{L, SXTNode{L, T}}}
end

function SXTree(n_nodes::L, ::Type{T} = Float64) where {L<:Integer, T}
    sxtree = SXTree(Vector{SXTNode{L, T}}(n_nodes), Dict{L, SXTNode{L, T}}[])

    for l in 1:n_nodes
        sxtree[l] = SXTNode(l, 0, sxtree)
    end

    sxtree
end

Base.show(io::IO, sxt::SXTree{L, T}) where {L, T} =
    print(io, "SXTree{$L, $T} with $(length(sxt.children)) vertices.")

# ============================================================================ #
# ITERATORS
# ============================================================================ #
mutable struct CousinIter{L, T}
    start::SXTNode{L, T}
    first::Bool
end

Base.start(ci::CousinIter) = ci.start
Base.done(ci::CousinIter, nd::SXTNode) = !ci.first && nd === ci.start
Base.next(ci::CousinIter, nd::SXTNode) =
    if ci.first
        ci.first = false
        nd, nd.next[]
    else
        nd, nd.next[]
    end

Base.iteratorsize(::Type{CousinIter}) = Base.SizeUnknown()
Base.iteratoreltype(::Type{CousinIter}) = Base.HasEltype()
Base.eltype(::Type{CousinIter{L, T}}) where {L, T} = SXTNode{L, T}

Base.show(io::IO, ci::CousinIter) =
    show(io, "CousinIter($(ci.start))")

"""
Iterate over all nodes with the same label on the same level.
"""
cousins(nd::SXTNode) = CousinIter(nd, true)

# ---------------------------------------------------------------------------- #
struct UpIter{L, T}
    start::SXTNode{L, T}
end
Base.start(ui::UpIter) = ui.start
Base.done(ui::UpIter, nd::SXTNode) = false
Base.done(ui::UpIter, sxt::SXTree) = true
Base.next(ui::UpIter, nd::SXTNode) = nd, nd.parent

Base.iteratorsize(::Type{UpIter}) = Base.HasLength()
Base.iteratoreltype(::Type{UpIter}) = Base.HasEltype()
Base.eltype(::Type{UpIter{L, T}}) where {L, T} = SXTNode{L, T}
Base.length(ui::UpIter) = ui.start.depth + 1

Base.show(io::IO, ci::UpIter) =
    show(io, "UpIter($(ci.start))")

"""
Iterate ancestors upwards until reaching the root node.
"""
ancestors(nd::SXTNode) = UpIter(nd)

# ---------------------------------------------------------------------------- #
struct BFIter{L, T}
    start::AbstractNode{L, T}
end

function Base.start(bi::BFIter{L, T}) where {L, T}
    q = Queue(SXTNode{L, T})
    foreach(n -> enqueue!(q, n), children(bi.start))
    q
end
Base.done(::BFIter, q) = isempty(q)
function Base.next(::BFIter, q)
    nd = dequeue!(q)
    foreach(n -> enqueue!(q, n), children(nd))
    nd, q
end
Base.show(bi::BFIter) = println("BFIter($(bi.start))")

Base.iteratorsize(::Type{BFIter}) = Base.SizeUnknown()
Base.iteratoreltype(::Type{BFIter}) = Base.HasEltype()
Base.eltype(::Type{BFIter{L, T}}) where {L, T} = SXTNode{L, T}

struct DFIter{L, T}
    start::AbstractNode{L, T}
end

function Base.start(di::DFIter{L, T}) where {L, T}
    s = Stack(SXTNode{L, T})
    foreach(n -> push!(s, n), (reverse ∘ collect ∘ children)(di.start))
    s
end
Base.done(::DFIter, s) = isempty(s)
function Base.next(::DFIter, s)
    nd = pop!(s)
    foreach(n -> push!(s, n), (reverse ∘ collect ∘ children)(nd))
    nd, s
end
Base.show(bi::DFIter) = println("DFIter($(bi.start))")

Base.iteratorsize(::Type{DFIter}) = Base.SizeUnknown()
Base.iteratoreltype(::Type{DFIter}) = Base.HasEltype()
Base.eltype(::Type{DFIter{L, T}}) where {L, T} = SXTNode{L, T}

"""
Traverse successors of node. Us BFIter for breadth first search (default) and
DFIter for depth first search.
"""
successors(an::AbstractNode, ::Type{BFIter}) = BFIter(an)
successors(an::AbstractNode) = BFIter(an)
successors(an::AbstractNode, ::Type{DFIter}) = DFIter(an)

"""
Iterate over children of a node.
"""
children(nd::SXTNode) = values(nd.children)
children(sxt::SXTree) = sxt.children

# ============================================================================ #
# GETTERS, SETTERS, UTILS
# ============================================================================ #
Base.getindex(an::AbstractNode, i) = an.children[i]
function Base.getindex(an::AbstractNode, is...)
    nd = an.children[is[1]]
    for i in is[2:end]
        nd = nd.children[i]
    end
    nd
end
Base.setindex!(nd::AbstractNode, m, i) = nd.children[i] = m

label(::SXTree{L, T}) where {L, T} = zero(L)
label(nd::SXTNode) = nd.label

data(nd::SXTNode) = get(nd.data)
n_vertices(sxt::SXTree) = length(sxt.children)

Base.vec(nd::SXTNode) = reverse!(collect(label(n) for n in ancestors(nd)))

dim(sxt::SXTree) = length(sxt.cousin_list)

function n_simplices(sxt::SXTree)
    n = 0
    for _ in successors(sxt)
        n += 1
    end
    n
end

# ============================================================================ #
# PRIVATE HELPERS
# ============================================================================ #
parent(nd::SXTNode) = nd.parent

has_child(sxt::SXTree, l) = l <= length(sxt.children)
has_child(nd::SXTNode, l) = haskey(nd.children, l)
has_children(an::AbstractNode, l) = all(has_child.(an, l))
num_children(an::AbstractNode) = length(an.children)
has_children(an::AbstractNode) = num_children(an) > 0

depth(nd::SXTNode) = nd.depth

function connect_cousins!(n1::SXTNode, n2::SXTNode)
    n1_next = n1.next[]
    n1.next[] = n2
    n1_next.prev[] = n2
    n2.prev[] = n1
    n2.next[] = n1_next
end

function add_child!(sxt::SXTree, nd::SXTNode, l)
    d = depth(nd) + 1
    new = SXTNode(l, d, nd)
    nd[l] = new

    # Connect to cousins or add a new entry in cousin list.
    if haskey(sxt.cousin_list[d], l)
        connect_cousins!(sxt.cousin_list[d][l], new)
    else
        sxt.cousin_list[d][l] = new
    end

    new
end

# ============================================================================ #
# OPERATIONS
# ============================================================================ #
"""
    insert_simplex!(sxtree::SXTree, simplex)

Insert a simplex into sxtree.
The simplex is given as an array of labels<:Integer.
"""
function insert_simplex!(sxtree::SXTree{L, T},
                         simplex::AbstractArray{L, 1}) where {L, T}

    simplex = sort(simplex)
    check_simplex(sxtree, simplex)

    while dim(sxtree) < length(simplex)-1
        push!(sxtree.cousin_list, Dict())
    end

    q = Queue(Tuple{SXTNode{L, T}, typeof(simplex)})
    for i in 2:length(simplex)
        enqueue!(q, (sxtree[simplex[i-1]], simplex[i:end]))
    end

    while !isempty(q)
        nd, sx = dequeue!(q)
        for i in 1:length(sx)
            if !has_child(nd, sx[i])
                add_child!(sxtree, nd, sx[i])
            end
            enqueue!(q, (nd[sx[i]], sx[i+1:end]))
        end
    end
end

"""
    is_face(sx1, sx2)

Check whether sx1 is a face of sx2.
"""
function is_face(sx::AbstractArray, nd::SXTNode)
    sx = sort(sx)
    check_simplex(sx)

    i = endof(sx)
    for p in ancestors(nd)
        i == 0 && break
        i > p.depth + 1 && break

        if label(p) == sx[i]
            i -= 1
        end
    end

    i == 0
end

is_face(nd1::SXTNode, nd2::SXTNode) = is_face(vec(nd1), nd2)
is_face(sx1::AbstractArray, sx2::AbstractArray) = issubset(sx1, sx2)
is_face(nd::SXTNode, sx::AbstractArray) = issubset(vec(nd), sx)


"""
    cofaces(sxtree, sx)

Find all cofaces of simplex `sx`, including itself.
"""
function cofaces(sxtree::SXTree{L, T}, sx::AbstractArray{L}) where {L, T}
    sx = sort(sx)
    check_simplex(sxtree, sx)

    d = length(sx) - 1
    cofaces = SXTNode{L, T}[]

    dopt = d
    for list in sxtree.cousin_list[d:end]
        if haskey(list, sx[end])
            for nd in cousins(list[sx[end]])
                if is_face(sx, nd)
                    push!(cofaces, nd)
                    foreach(n -> push!(cofaces, n), successors(nd))
                end
            end
        end
    end

    cofaces
end

"""
Remove node WITHOUT fixing cofaces or SXTree dimension.
Do not use this, use rem_simplex! instead.
"""
function rem_node!(sxtree::SXTree{L, T}, nd::SXTNode{L, T}) where {L, T}
    depth(nd) == 0 && throw(ArgumentError("Can't remove top nodes!"))

    cousins = sxtree.cousin_list[depth(nd)]
    if nd.next[] === nd
        # Only cousin.
        delete!(sxtree.cousin_list[depth(nd)], label(nd))
    else
        # Remove from list.
        cousins[label(nd)] === nd && (cousins[label(nd)] = nd.next[])
        nd.next[].prev[] = nd.prev[]
        nd.prev[].next[] = nd.next[]
    end

    delete!(parent(nd).children, label(nd))
end

"""
    rem_simplex!(sxtree, sx)

Remove simplex `sx` from `sxtree`. `sx` can be given as a node or as an array.
"""
function rem_simplex!(sxtree::SXTree{L, T}, sx::AbstractArray{L}) where {L, T}
    sort(sx)
    check_simplex(sxtree, sx)

    foreach(n -> rem_node!(sxtree, n), cofaces(sxtree, sx))

    d = dim(sxtree)
    while isempty(sxtree.cousin_list[d])
        pop!(sxtree.cousin_list)
        d -= 1
    end
end

function rem_simplex!(sxtree::SXTree{L, T}, nd::SXTNode{L, T}) where {L, T}
    rem_simplex!(sxtree, vec(nd))
end
