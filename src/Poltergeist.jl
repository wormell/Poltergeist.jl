module Poltergeist
  using Base, ApproxFun, BandedMatrices, Compat, DualNumbers, StaticArrays, IntervalSets#, PyPlot#, FastTransforms

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,
>=,./,/,.^,^,\,∪,transpose, size, length,  eltype, inv, mod, convert#issymmetric,
import ApproxFun: domainspace, rangespace, domain, israggedbelow, RaggedMatrix, resizedata!, colstop, CachedOperator, Infinity,
fromcanonicalD, tocanonicalD, fromcanonical, tocanonical, space

export (..)

## TEMPORARY PENDING APPROXFUN UPDATE
temp_in(x,dom) = (x >= dom.a) && (x <= dom.b)

include("general.jl")
include("MarkovBranch.jl")
include("MarkovMap.jl")
include("CircleMap.jl")
include("Transfer.jl")
include("InducedTransfer.jl")
include("Solution.jl")
include("TimeSeries.jl")
include("poetry.jl")

end # module
