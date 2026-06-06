module Ellasy

# Write your package code here.

import LinearAlgebra

include("./matrices.jl")
include("./ladder.jl")

export SingleLadder
export getLadderType
export getCircuitTypes, buildSecondOrderSystem

# New naming
export secondOrder, firstOrder, firstOrderSubmatrices

# Deprecated naming
export buildFOSubmatrices, buildFirstOrderSystem

end