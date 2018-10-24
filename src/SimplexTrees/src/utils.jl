function check_simplex(sx)
    if !(eltype(sx) <: Integer)
        throw(ArgumentError("Invalid simplex '$sx'!" *
                            "Labels should <: Integer!"))
    end
    if isa(sx, AbstractSet)
        sort!(collect(sx))
    else
        if !allunique(sx)
            throw(ArgumentError("Invalid simplex '$sx'!" *
                                "Vertices in simplex should be unique."))
        end
        if !issorted(sx)
            throw(ArgumentError("Invalid simplex '$sx'!" *
                                "The labels in simplex should be sorted."))
        end
        sx
    end
end

function check_simplex(sxtree, sx)
    check_simplex(sx)
    if any(sx .> num_children(sxtree))
        lbl = filter(l -> l > num_children(sxtree), sx)
        throw(ArgumentError("Invalid simplex '$sx'!" *
                            "Label $lbl too large for this tree."))
    end
end
