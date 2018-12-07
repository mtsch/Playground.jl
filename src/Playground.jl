module Playground

using Reexport
include("geodesiccomplex.jl")
@reexport using Playground.GeodesicComplexes
include("equilateraliterator.jl")
@reexport using Playground.EquilateralIterators
include("persistence.jl")

export persistence, cleanup

end
