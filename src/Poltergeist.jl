__precompile__()

module Poltergeist
<<<<<<< HEAD
  using Base, ApproxFun, BandedMatrices, Compat, DualNumbers, StaticArrays, IntervalSets#, PyPlot#, FastTransforms
  # VERSION < v"0.7" ? (using Compat.LinearAlgebra) : (using LinearAlgebra)
  using LinearAlgebra

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,
  >=,./,/,.^,^,\,∪,∘,transpose, size, length,  eltype, inv, mod, convert#issymmetric,
@compat import LinearAlgebra: eigvals, eigvecs
import ApproxFun: domainspace, rangespace, domain, israggedbelow, RaggedMatrix,
  resizedata!, colstop, CachedOperator, Infinity, IntervalDomain,
  fromcanonicalD, tocanonicalD, fromcanonical, tocanonical, space, eigs, cfstype
=======
  using Base, ApproxFun, BandedMatrices, Compat, DualNumbers
  using StaticArrays, IntervalSets, LightGraphs#, PyPlot#, FastTransforms

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,
  >=,./,/,.^,^,\,∪,∘,transpose, size, length,  eltype, inv, mod, convert,#issymmetric,
  eigvals, eigs, show
import ApproxFun: domainspace, rangespace, domain, israggedbelow, RaggedMatrix,
  resizedata!, colstop, CachedOperator, Infinity, IntervalDomain, UnionDomain,
  fromcanonicalD, tocanonicalD, fromcanonical, tocanonical, space, prectype
>>>>>>> 4c33d1554aa02b6a20399f0090f695de8ad3e517

export (..), Segment, PeriodicInterval, Fun

## TEMPORARY PENDING APPROXFUN UPDATE
temp_in(x,dom) = (x >= dom.a) && (x <= dom.b)

include("general.jl")
include("ExpandingBranch.jl")
include("IntervalMap.jl")
include("CircleMap.jl")
include("ComposedMaps.jl")
include("HGraphs.jl")
include("hofbauerextension.jl")
include("InducedMap.jl")
include("BackSum.jl")
include("Transfer.jl")
include("Solution.jl")
include("TimeSeries.jl")
include("conversion.jl")
include("conjugacy.jl")
include("statistics.jl")
include("poetry.jl")
include("examples.jl")

end # module
