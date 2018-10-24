# ============================================================================ #
# DictSet
# ============================================================================ #
"""
Basic dict-like thing that holds multiple values for each key in a set.
"""
struct DictSet{K, V}
    data::Dict{K, Set{V}}
end

DictSet{K, V}() where {K, V} = DictSet{K, V}(Dict{K, Set{V}}())

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

function Base.delete!(ds::DictSet{K}, k::K) where K
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

function Base.delete!(ds::DictSet{K, V}, vals::AbstractSet{V}) where {K, V}
    for (k, v) in ds.data
        setdiff!(v, vals)
        if isempty(v)
            delete!(ds.data, k)
        end
    end
    ds
end

Base.keys(ds::DictSet) = keys(ds.data)
Base.values(ds::DictSet) = union(values(ds.data)...)
Base.haskey(ds::DictSet) = haskey(ds.data)
Base.isempty(ds::DictSet) = isempty(ds.data)

# ============================================================================ #
"""
  validatesimplex(σ)

Check if a simplex σ is valid. Throw an ArgumentError if not.
A simplex should contain Integers, be sorted and contain no duplicates.
"""
function validatesx(σ)
    if !(eltype(σ) <: Integer)
        throw(ArgumentError("Invalid simplex '$σ'!" *
                            "Labels should be Integers!"))
    end
    if !allunique(σ)
        throw(ArgumentError("Invalid simplex '$σ'!" *
                            "Vertices in simplex should be unique."))
    end
    if !issorted(σ)
        throw(ArgumentError("Invalid simplex '$σ'!" *
                            "The labels in simplex should be sorted."))
    end
    σ
end

"""
  isface(τ, σ)

Check if τ is a face of σ. Assumes the simplices are valid.
"""
function isface(τ, σ)
    i = 1; j = 1
    while i <= length(τ) && j <= length(σ)
        if τ[i] > σ[j]
            j += 1
        elseif τ[i] == σ[j]
            i += 1; j += 1
        else
            return false
        end
    end
    i == length(τ) + 1
end
