module EquilateralIterators
export equilaterals, getsortededges

"""
it iterates equilateral triangles!
"""
struct EquilateralIterator{T, M<:AbstractMatrix{T}}
    dists       ::M
    sortededges ::Vector{Tuple{Int, Int}}
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
equilaterals(dists::AbstractMatrix{T}, tol) where T =
    EquilateralIterator{T, typeof(dists)}(dists, getsortededges(dists), T(tol))

"""
    getsortededges(dists)

Get ordering on edges (values in matrix).
"""
getsortededges(dists) =
    sort!([(i, j) for i in 2:size(dists, 1) for j in 1:i-1],
          lt = (e1, e2) -> isless(dists[e1...], dists[e2...]))

function Base.iterate(ei::EquilateralIterator, st = (1, 1))
    dists = ei.dists
    sortededges = ei.sortededges
    tol = ei.tol

    eindex, kstart = st
    if kstart > size(dists, 1)
        kstart = 1
        eindex += 1
    end
    n = size(dists, 1)

    @inbounds while eindex ≤ length(sortededges)
        i, j = sortededges[eindex]
        a = dists[i, j]

        for k in kstart:n
            (k == i || k == j) && continue
            b = dists[k, i]
            c = dists[k, j]
            # avoid problems when points are not in general position.
            a > b || k > i || continue
            a > c || k > j || continue

            if 0 ≤ a - b ≤ tol && 0 ≤ a - c ≤ tol && abs(b - c) ≤ tol
                return (i, j, k), (eindex, k + 1)
            end
        end
        kstart = 1
        eindex += 1
    end

    nothing
end

end
