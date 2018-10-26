using LinearAlgebra

struct EquilateralIterator{T, M<:AbstractMatrix{T}}
    dists       ::M
    sortededges ::Vector{Tuple{Int, Int}}
    indexmatrix ::Symmetric{Int, Matrix{Int}} # <- remove
    tol         ::T
end

Base.IteratorSize(::EquilateralIterator) = Base.SizeUnknown()
Base.IteratorEltype(::EquilateralIterator) = Base.HasEltype()
Base.eltype(::EquilateralIterator) = Tuple{Int, Int, Int}
Base.eltype(::Type{EquilateralIterator}) = Tuple{Int, Int, Int}

"""
    equilaterals(dists, tol)

Iterator that returns approximately (up to `tol`) equilateral triangles in
distance matrix `dists`, sorted by longest side length.
"""
function equilaterals(dists::AbstractMatrix{T}, tol) where T
    t = T(tol)
    sortededges, indexmatrix = getedgeorder(M)
    EquilateralIterator{T, typeof(M)}(dists, sortededges, indexmatrix, t)
end

"""
    getedgeorder(dists)

Get ordering on edges (values in matrix). Returns sorted edges and index matrix.
"""
function getedgeorder(dists)
    n = size(dists, 1)
    sortededges = sort!([(i, j) for i in 2:n for j in 1:i-1],
                        lt = (e1, e2) -> isless(dists[e1...], dists[e2...]))
    indexmatrix = fill(0, n, n)

    for (i, e) in enumerate(sortededges)
        indexmatrix[e...] = i
    end
    sortededges, Symmetric(indexmatrix, :L)
end

function Base.iterate(ei::EquilateralIterator, st = (1, 1))
    dists = ei.dists
    sortededges = ei.sortededges
    indexmatrix = ei.indexmatrix
    tol = ei.tol

    estart, kstart = st
    n = size(dists, 1)

    @inbounds for eindex in estart:length(sortededges)
        i, j = sortededges[eindex]
        a = dists[i, j]

        for k in kstart:n
            (k == i || k == j) && continue
            b = M[k, i]
            c = M[k, j]

            if 0 < a - b ≤ tol && 0 < a - c ≤ tol && abs(b - c) ≤ tol
                return (i, j, k), (eindex, k + 1)
            end
        end
    end

    nothing
end
