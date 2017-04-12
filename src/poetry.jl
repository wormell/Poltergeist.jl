export acim, linearresponse, correlationsum, birkhoff_cov, birkhoff_var

acim(S::SolutionInvWrapper) = S.op\S.u
acim(L::Operator,u::Fun=uniform(domainspace(L))) = SolutionInv(L,u)\u
acim(M::AbstractMarkovMap) = acim(Transfer(M))

function normalise_for_correlation(f::Fun,r::Fun=uniform(space(f)))
  rf = r*f
  rf -= r*sum(rf)/sum(r)
end

linearresponse(S::SolutionInvWrapper,X::Fun) = S\(-acim(S)*X)'
function correlationsum(S::SolutionInvWrapper,A::Fun,r=acim(S))
  S \ normalise_for_correlation(A,r)
end

function birkhoffcov(S::SolutionInvWrapper,A::Fun,B::Fun)
  r = acim(S)
  sum(B*correlationsum(S,A,r)) + sum(A*correlationsum(S,B,r)) - sum(B*normalise_for_correlation(A,r))
end

function birkhoffvar(S::SolutionInvWrapper,A::Fun)
  r = acim(S)
  2sum(A*correlationsum(S,A,r)) - sum(A*normalise_for_correlation(A,r))
end

for OP in (:linearresponse,:correlationsum,:birkhoffvar)
  @eval $OP(L::Operator,X::Fun) = $OP(SolutionInv(L),X)
  @eval $OP(M::AbstractMarkovMap,X::Fun) = $OP(Transfer(M),X)
end
birkhoffcov(L::Operator,X::Fun,Y::Fun) = birkhoffcov(SolutionInv(L),X,Y)
birkhoffcov(M::AbstractMarkovMap,X::Fun,Y::Fun) = birkhoffcov(Transfer(M),X,Y)


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

immutable MarkovMapDerivative{M<:AbstractMarkovMap}
  m::M
end
@compat (md::MarkovMapDerivative)(x::Number) = mapD(md.m,x)
Base.ctranspose(b::AbstractMarkovMap) = MarkovMapDerivative(b)


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
