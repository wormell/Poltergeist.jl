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

acim(S::SchurInvWrapper) = S.op\S.u
acim(L::Operator,u::Fun=uniform(domainspace(L))) = SchurInv(L,u)\u
acim(M::MarkovMap) = acim(Transfer(M))

linearresponse(S::SchurInvWrapper,X::Fun) = S\(-acim(S)*X)'
correlationsum(S::SchurInvWrapper,A::Fun) = S\(acim(S)*A)

for OP in (:linearresponse,:correlationsum)
  @eval $OP(L::Operator,X::Fun) = $OP(SchurInv(L),X)
  @eval $OP(M::MarkovMap,X::Fun) = $OP(Transfer(M),X)
end

# Markov derivatives
immutable MarkovBranchDerivative{B<:MarkovBranch}
  b::B
end
(bd::MarkovBranchDerivative)(x::Number) = mapD(bd.b,x)
ctranspose(b::MarkovBranch) = MarkovBranchDerivative(b)

immutable MarkovBranchInverse{B<:MarkovBranch}
  b::B
end
(bi::MarkovBranchInverse)(x::Number) = mapinv(bd.b,x)
inv(b::MarkovBranch) = MarkovBranchInverse(b)
## TODO one day: convert routine to MarkovBranch?

immutable MarkovBranchDerivativeInverse{B<:MarkovBranch}
  b::B
end
(bi::MarkovBranchDerivativeInverse)(x::Number) = mapinvD(bd.b,x)
ctranspose(bi::MarkovBranchInverse) = MarkovBranchDerivativeInverse(bi.b)

