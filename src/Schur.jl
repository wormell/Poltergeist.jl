export SchurInv

# Schur wrappers
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
SchurInv(M::MarkovMap,u::Fun=uniform(Space(domain(M)))) = SchurInv(Transfer(M,space(u)),u)
