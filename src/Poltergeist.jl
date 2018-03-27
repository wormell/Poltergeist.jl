module Poltergeist
  using Base, ApproxFun, BandedMatrices, Compat, DualNumbers, StaticArrays, IntervalSets#, PyPlot#, FastTransforms

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,
  >=,./,/,.^,^,\,∪,∘,transpose, size, length,  eltype, inv, mod, convert#issymmetric,
import ApproxFun: domainspace, rangespace, domain, israggedbelow, RaggedMatrix,
  resizedata!, colstop, CachedOperator, Infinity, IntervalDomain,
  fromcanonicalD, tocanonicalD, fromcanonical, tocanonical, space

export (..), Segment, PeriodicInterval, Fun

## TEMPORARY PENDING APPROXFUN UPDATE
temp_in(x,dom) = (x >= dom.a) && (x <= dom.b)

include("general.jl")
include("ExpandingBranch.jl")
include("MarkovMap.jl")
include("CircleMap.jl")
include("ComposedMaps.jl")
include("Transfer.jl")
include("Solution.jl")
include("TimeSeries.jl")
include("conjugacy.jl")
include("statistics.jl")
include("poetry.jl")
include("examples.jl")

end # module
