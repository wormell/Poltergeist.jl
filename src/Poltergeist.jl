__precompile__()

module Poltergeist

using Base, ApproxFun, BandedMatrices, Compat, DualNumbers, StaticArrays, IntervalSets, LightGraphs#, PyPlot#, FastTransforms
using LinearAlgebra

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,
  >=,./,/,.^,^,\,∪,∘,transpose, size, length,  eltype, inv, mod, convert#issymmetric,
@compat import LinearAlgebra: eigvals, eigvecs
import ApproxFun: domainspace, rangespace, domain, israggedbelow, RaggedMatrix,
  resizedata!, colstop, CachedOperator, Infinity, IntervalDomain, UnionDomain,
  fromcanonicalD, tocanonicalD, fromcanonical, tocanonical, space, eigs, cfstype, prectype

export (..), Segment, PeriodicInterval, Fun, eigs

## TEMPORARY PENDING APPROXFUN UPDATE
temp_in(x,dom) = (x >= dom.a) && (x <= dom.b)

include("general.jl")
include("maps/maps.jl")
include("HGraphs.jl")
include("hofbauerextension.jl")
include("InducedMap.jl")
include("BackSum.jl")
include("spectral/spectral.jl")
include("timeseries.jl")
include("poetry.jl")

end # module
