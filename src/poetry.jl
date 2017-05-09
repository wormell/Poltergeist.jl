export acim, linearresponse, correlationsum, birkhoffcov, birkhoffvar

acim(S::SolutionInvWrapper) = S.op\S.u
acim(L::Operator,u::Fun=uniform(domainspace(L))) = SolutionInv(L,u)\u
acim(M::AbstractMarkovMap) = acim(Transfer(M))

function zero_to(f::Fun,r::Fun=uniform(space(f)))
  rf = r*f
  rf -= r*sum(rf)/sum(r)
end

function linearresponse(S::SolutionInvWrapper,X::Fun)
  @assert all(isapprox.(X.(∂(domain(X))),0))
  S\(-acim(S)*X)'
end
correlationsum(S::SolutionInvWrapper,A::Fun,r=acim(S)) = S \ zero_to(A,r)

function birkhoffcov(S::SolutionInvWrapper,A::Fun,B::Fun)
  r = acim(S); Az = zero_to(A,r); Bz = zero_to(B,r)
  sum(B*(S\Az)) + sum(A*(S\Bz)) - sum(B*Az)
end

function birkhoffvar(S::SolutionInvWrapper,A::Fun)
  Az = zero_to(A,acim(S))
  2sum(A*(S\Az)) - sum(A*Az)
end

#TODO: Lyapunov exponent
# function lyapunov(S::SolutionInvWrapper{})

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
