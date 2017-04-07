export acim, linearresponse, correlationsum, birkhoff_cov, birkhoff_var

acim(S::SchurInvWrapper) = S.op\S.u
acim(L::Operator,u::Fun=uniform(domainspace(L))) = SchurInv(L,u)\u
acim(M::AbstractMarkovMap) = acim(Transfer(M))

function normalise_for_correlation(f::Fun,r::Fun=uniform(space(f)))
  rf = r*f
  rf -= r*sum(rf)/sum(r)
end

linearresponse(S::SchurInvWrapper,X::Fun) = S\(-acim(S)*X)'
function correlationsum(S::SchurInvWrapper,A::Fun)
  r = acim(S)
  S \ normalise_for_correlation(A,r)
end

birkhoff_cov(S::SchurInvWrapper,A::Fun,B::Fun) =
  sum(B*correlationsum(S,A)) + sum(A*correlationsum(S,B)) - sum(B*normalise_for_correlation(A,acim(S)))
birkhoff_var(S::SchurInvWrapper,A::Fun) =
  2sum(A*correlationsum(S,A)) - sum(A*normalise_for_correlation(A,acim(S)))

for OP in (:linearresponse,:correlationsum,:birkhoff_var)
  @eval $OP(L::Operator,X::Fun) = $OP(SchurInv(L),X)
  @eval $OP(M::AbstractMarkovMap,X::Fun) = $OP(Transfer(M),X)
end
birkhoff_cov(L::Operator,X::Fun,Y::Fun) = birkhoff_cov(SchurInv(L),X,Y)
birkhoff_cov(M::AbstractMarkovMap,X::Fun,Y::Fun) = birkhoff_cov(Transfer(M),X,Y)


# Markov derivatives
immutable MarkovBranchDerivative{B<:MarkovBranch}
  b::B
end
@compat (bd::MarkovBranchDerivative)(x::Number) = mapD(bd.b,x)
Base.ctranspose(b::MarkovBranch) = MarkovBranchDerivative(b)
# Base.transpose(b::MarkovBranch) = MarkovBranchDerivative(b)

immutable MarkovBranchInverse{B<:MarkovBranch}
  b::B
end
@compat (bi::MarkovBranchInverse)(x::Number) = mapinv(bi.b,x)
inv(b::MarkovBranch) = MarkovBranchInverse(b)
## TODO one day: convert routine to MarkovBranch?

immutable MarkovBranchDerivativeInverse{B<:MarkovBranch}
  b::B
end
@compat (bi::MarkovBranchDerivativeInverse)(x::Number) = mapinvD(bd.b,x)
Base.ctranspose(bi::MarkovBranchInverse) = MarkovBranchDerivativeInverse(bi.b)
# Base.transpose(b::MarkovBranch) = MarkovBranchDerivative(b)

immutable MarkovMapDerivative{M<:MarkovMap}
  m::M
end
@compat (md::MarkovMapDerivative)(x::Number) = mapD(md.m,x)
Base.ctranspose(b::MarkovMap) = MarkovMapDerivative(b)


# # plotting
# function plot(m::MarkovMap)
#   pts = eltype(m)[]
#   vals = eltype(m)[]
#   sp = sortperm([minimum(∂(domain(b))) for b in branches(m)])
#   for b in branches(m)[sp]
#     append!(pts,points(domain(b),100))
#     append!(vals,broadcast(b,points(domain(b),100)))
#   end
#   p = PyPlot.plot(pts,vals)
#   xlim(∂(domain(m)))
#   ylim(∂(rangedomain(m)))
#   xlabel("\$x\$")
#   ylabel("\$f(x)\$")
#   p
# end

# function plot(m::AbstractMarkovMap)
#   pts = points(domain(m),500)
#   vals = broadcast(m,pts)
#   p = PyPlot.plot(pts,vals)
#   xlim(∂(domain(m)))
#   ylim(∂(rangedomain(m)))
#   xlabel("\$x\$")
#   ylabel("\$f(x)\$")
#   p
# end
