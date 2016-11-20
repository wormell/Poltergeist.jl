export acim, linearresponse, correlationsum

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
Base.ctranspose(b::MarkovBranch) = MarkovBranchDerivative(b)

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
Base.ctranspose(bi::MarkovBranchInverse) = MarkovBranchDerivativeInverse(bi.b)

