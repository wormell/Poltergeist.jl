module Poltergeist
  using Base, ApproxFun, BandedMatrices, Compat, DualNumbers, FixedSizeArrays#, PyPlot#, FastTransforms

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,
>=,./,/,.^,^,\,∪,transpose, size, length,  eltype, inv#issymmetric,
import ApproxFun: domainspace, rangespace, domain, israggedbelow, RaggedMatrix, resizedata!, colstop, CachedOperator, Infinity,
fromcanonicalD, tocanonicalD, default_raggedmatrix, fromcanonical, tocanonical, space

include("general.jl")
include("AbelFunction.jl")
include("MarkovMap.jl")
include("Transfer.jl")
include("InducedTransfer.jl")
include("Schur.jl")
include("TimeSeries.jl")
include("poetry.jl")
end # module
