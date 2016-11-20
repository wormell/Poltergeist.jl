export acim, linearresponse, correlationsum

acim(S::SchurInvWrapper) = S.op\S.u
acim(L::Operator,u::Fun=uniform(domainspace(L))) = SchurInv(L,u)\u
acim(M::AbstractMarkovMap) = acim(Transfer(M))

linearresponse(S::SchurInvWrapper,X::Fun) = S\(-acim(S)*X)'
correlationsum(S::SchurInvWrapper,A::Fun) = S\(acim(S)*A)

for OP in (:linearresponse,:correlationsum)
  @eval $OP(L::Operator,X::Fun) = $OP(SchurInv(L),X)
  @eval $OP(M::AbstractMarkovMap,X::Fun) = $OP(Transfer(M),X)
end

# Markov derivatives
immutable MarkovBranchDerivative{B<:MarkovBranch}
  b::B
end
(bd::MarkovBranchDerivative)(x::Number) = mapD(bd.b,x)
Base.ctranspose(b::MarkovBranch) = MarkovBranchDerivative(b)
# Base.transpose(b::MarkovBranch) = MarkovBranchDerivative(b)

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
# Base.transpose(b::MarkovBranch) = MarkovBranchDerivative(b)

immutable MarkovMapDerivative{M<:MarkovMap}
  m::M
end
(md::MarkovMapDerivative)(x::Number) = mapD(md.m,x)
Base.ctranspose(b::MarkovMap) = MarkovMapDerivative(b)


# plotting
function plot(m::MarkovMap)
  pts = eltype(m)[]
  vals = eltype(m)[]
  sp = sortperm([minimum(∂(domain(b))) for b in branches(m)])
  for b in branches(m)[sp]
    append!(pts,points(domain(b),100))
    append!(vals,broadcast(b,points(domain(b),100)))
  end
  p = PyPlot.plot(pts,vals)
  xlim(∂(domain(m)))
  ylim(∂(rangedomain(m)))
  xlabel("\$x\$")
  ylabel("\$f(x)\$")
  p
end

function plot(m::AbstractMarkovMap)
  pts = points(domain(m),500)
  vals = broadcast(m,pts)
  p = PyPlot.plot(pts,vals)
  xlim(∂(domain(m)))
  ylim(∂(rangedomain(m)))
  xlabel("\$x\$")
  ylabel("\$f(x)\$")
  p
end
