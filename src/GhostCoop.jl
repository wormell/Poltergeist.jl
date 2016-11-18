module GhostCoop
  using Base, ApproxFun, BandedMatrices, Compat#, FastTransforms

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,
>=,./,/,.^,^,\,∪,transpose, size, length, issymmetric, eltype#, maximum, minimum
import ApproxFun: domainspace, rangespace, domain, israggedbelow, RaggedMatrix, resizedata!, colstop, CachedOperator, Infinity,
fromcanonicalD, tocanonicalD, default_raggedmatrix

include("general.jl")
include("MarkovMap.jl")
include("Transfer.jl")

type SchurInvWrapper{QR<:ApproxFun.QROperator,F<:Fun,T} <: Operator{T}
  op::QR
  u::F
end
SchurInvWrapper(op::ApproxFun.QROperator,f::Fun) = SchurInvWrapper{typeof(op),typeof(f),eltype(op)}(op,f)
ApproxFun.qrfact!(op::SchurInvWrapper) = op
ApproxFun.qrfact(op::SchurInvWrapper) = op
#SchurInvWrapper(op::Operator) = SchurInvWrapper(op,uniform(domainspace(op)))
ApproxFun.@wrapper SchurInvWrapper
ApproxFun.linsolve(L::SchurInvWrapper,b;kwds...) = linsolve(L.op,b;kwds...)

uniform(S::Space) = Fun(1.,S)/sum(Fun(1.,S))
uniform(D::Domain) = Fun(1.,D)/sum(Fun(1.,D))
function SchurInv(L::Operator,u::Fun=uniform(domainspace(L)))
  @assert domain(L) == rangedomain(L)
  SchurInvWrapper(qrfact(I-L + cache(u/sum(u) * DefiniteIntegral(domainspace(L)))),u)
end

acim(L::Operator,u::Fun=uniform(domainspace(L))) = SchurInv(L,u)\u
acim(S::SchurInvWrapper) = S.op\S.u
linearresponse(S::SchurInvWrapper,X::Fun) = S\(-acim(S)*X)'
linearresponse(L::Operator,X::Fun) = linearresponse(SchurInv(L),X)
correlationsum(S::SchurInvWrapper,A::Fun) = S\(acim(S)*A)
correlationsum(L::Operator,A::Fun) = correlationsum(SchurInv(L),A)
end # module
