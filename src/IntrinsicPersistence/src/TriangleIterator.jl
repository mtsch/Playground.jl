struct Triangle{T}
    vertices ::NTuple{3, Int}
    r        ::T
end

Triangle(vs, r::T) where T = Triangle{T}(vs, r)

function Base.show(io::IO, Δ::Triangle{T}) where T
    if T <: Number
        print(io, "Δ($(Δ.vertices), $(round(Δ.r, 5)))")
    else
        print(io, "Δ($(Δ.vertices), $(Δ.r))")
    end
end

# Return "empty" triangle. Used as a placeholder value when the iterator is done.
nulltriangle(::Type{T}) where T = Triangle{T}((0, 0, 0), zero(T))

# ============================================================================ #
# Equilateral Iterator
# ============================================================================ #

# Return a vector that maps ordering to indices in distance matrix, for example,
# ordertoedge(M)[3] returns the third shortest edge (=third smallest entry) in M
function ordertoedge(M)
    n = size(M, 1)
    res = Vector{Tuple{Int, Int}}(n * (n - 1) ÷ 2)
    i = 1
    for j in 2:n, k in 1:j-1
        res[i] = (j, k)
        i += 1
    end
    sort!(res, lt = (e1, e2) -> isless(M[e1...], M[e2...]))
    res
end

struct EquilateralIterator{T, M<:AbstractMatrix{T}}
    M     ::M
    edges ::Vector{Tuple{Int, Int}}
    tol   ::T
    upper ::T
end

"""
    equilaterals(M; tol = sampledensity(M), upper = Inf)

Iterator that returns approximately (up to `tol`) equilateral triangles in
distance matrix `M`, sorted by longest side length. Skips triangles with
longest side length above `upper`.
"""
function equilaterals(M::AbstractMatrix{T};
                      tol   = sampledensity(M),
                      upper = typemax(T)) where T
    t = convert(T, tol)
    u = convert(T, upper)
    edges = ordertoedge(M)
    EquilateralIterator{T, typeof(M)}(M, edges, t, u)
end

function Base.start(ei::EquilateralIterator{T}) where T
    M     = ei.M
    edges = ei.edges
    tol   = ei.tol
    upper = ei.upper

    @inbounds for eindex in eachindex(edges)
        i, j = edges[eindex]
        a = M[i, j]
        a > upper && return (true, nulltriangle(T), 0, 0)
        a ≤ 0 && continue

        for k in 1:size(M, 1)
            (k == i || k == j || i == j) && continue
            b = M[k, i]
            b ≤ 0 && continue
            c = M[k, j]
            c ≤ 0 && continue
            if 0 ≤ a - b ≤ tol && 0 ≤ a - c ≤ tol && abs(b - c) ≤ tol
                return (false, Triangle((i, j, k), a), eindex, k + 1)
            end
        end
    end
    (true, nulltriangle(T), 0, 0)
end

function Base.next(ei::EquilateralIterator{T}, st) where T
    M     = ei.M
    edges = ei.edges
    tol   = ei.tol
    upper = ei.upper

    _, Δ, eindex, kstart = st
    n = size(M, 1)

    if kstart > n
        eindex += 1
        kstart  = 1
    end

    if eindex > length(edges) || M[edges[eindex]...] > upper
        return Δ, (true, nulltriangle(T), 0, 0)
    end

    @inbounds while eindex ≤ length(edges)
        i, j = edges[eindex]
        a = M[i, j]
        a > upper && return Δ, (true, nulltriangle(T), 0, 0)
        a ≤ 0 && continue

        for k in kstart:n
            (k == i || k == j || i == j) && continue
            b = M[k, i]
            b ≤ 0 && continue
            c = M[k, j]
            c ≤ 0 && continue
            if 0 ≤ a - b ≤ tol && 0 ≤ a - c ≤ tol && abs(b - c) ≤ tol
                return Δ, (false, Triangle((i, j, k), a), eindex, k + 1)
            end
        end

        eindex += 1
        kstart  = 1
    end

    Δ, (true, nulltriangle(T), 0, 0)
end

Base.done(ei::EquilateralIterator, st) = st[1]

Base.iteratorsize(::EquilateralIterator)   = Base.SizeUnknown()
Base.iteratoreltype(::EquilateralIterator) = Base.HasEltype()

Base.eltype(::EquilateralIterator{T}) where T       = Triangle{T}
Base.eltype(::Type{EquilateralIterator{T}}) where T = Triangle{T}
