module Playground

using Reexport
include("geodesiccomplex.jl")
@reexport using Playground.GeodesicComplexes
include("triangleiterator.jl")
@reexport using Playground.EquilateralIterators

end
