module SimplexTrees

using DataStructures

export
    # Types/SXTree
    AbstractNode, SXTNode, SXTree, CousinIterator,

    # Types/SAL
    SALNode, SAL,

    # Iterators
    children, cousins, ancestors, successors, BFIter, DFIter,

    # Getters, setters
    dim, n_vertices, label, data, n_simplices,

    # Operations
    insert_simplex!, is_face, cofaces, rem_simplex!,

    isface

include("utils.jl")
include("simplex_tree.jl")
include("simplex_array_list2.jl")

end
