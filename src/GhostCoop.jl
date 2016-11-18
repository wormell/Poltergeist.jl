module GhostCoop
  using Base, ApproxFun, BandedMatrices, Compat#, FastTransforms

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,
>=,./,/,.^,^,\,∪,transpose, size, length, issymmetric, eltype#, maximum, minimum
import ApproxFun: domainspace, rangespace, domain, israggedbelow, RaggedMatrix, resizedata!, colstop, CachedOperator, Infinity,
fromcanonicalD, tocanonicalD, default_raggedmatrix

include("general.jl")
include("MarkovMap.jl")
include("Transfer.jl")
include("poetry.jl")
end # module
