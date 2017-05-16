export acim, linearresponse, correlationsum, birkhoffcov, birkhoffvar
export lyapunov
export ∘

(∘)(f,g) = x->f(g(x))

# Lanford map
for (fn,dfn,n) in ((:lanford_v0,:lanford_dv0,25),(:lanford_v1,:lanford_dv1,17))
  @eval ($fn)(x) = (5-sqrt($n-8x))/2
  @eval ($dfn)(x) = 2/sqrt($n-8x)
end
lanford_brk = lanford_v0(1)
lanford = MarkovMap([lanford_v0,lanford_v1],[0..lanford_brk,lanford_brk..1];dir=Reverse,diff=[lanford_dv0,lanford_dv1])

# Overloads

acim(K::SolutionInvWrapper) = K.op\K.u
acim(L::Operator,u::Fun=uniform(domainspace(L))) = SolutionInv(L,u)\u
acim(M::AbstractMarkovMap) = acim(Transfer(M))

function zero_to(f::Fun,r::Fun=uniform(space(f)))
  rf = r*f
  rf -= r*sum(rf)/sum(r)
end

function linearresponse(K::SolutionInvWrapper,X::Fun)
  @assert all(isapprox.(X.(∂(domain(X))),0;atol=sqrt(eps(eltype(X))*arclength(domain(X)))))
  K\(-acim(K)*X)'
end
correlationsum(K::SolutionInvWrapper,A::Fun,r=acim(K)) = K \ zero_to(A,r)

function birkhoffcov(K::SolutionInvWrapper,A::Fun,B::Fun)
  r = acim(K); Az = zero_to(A,r); Bz = zero_to(B,r)
  sum(B*(K\Az)) + sum(A*(K\Bz)) - sum(B*Az)
end

function birkhoffvar(K::SolutionInvWrapper,A::Fun)
  Az = zero_to(A,acim(K))
  2sum(A*(K\Az)) - sum(A*Az)
end

for OP in (:linearresponse,:correlationsum,:birkhoffvar)
  @eval $OP(L::Operator,X::Fun) = $OP(SolutionInv(L),X)
  @eval $OP(M::AbstractMarkovMap,X::Fun) = $OP(Transfer(M),X)
  @eval $OP(obj::Operator,X) = $OP(obj,Fun(X,domainspace(obj)))
  @eval $OP(obj::AbstractMarkovMap,X) = $OP(obj,Fun(X,domain(obj)))
end
birkhoffcov(L::Operator,X::Fun,Y::Fun) = birkhoffcov(SolutionInv(L),X,Y)
birkhoffcov(M::AbstractMarkovMap,X::Fun,Y::Fun) = birkhoffcov(Transfer(M),X,Y)

#Lyapunov exponent
#TODO: integrate into transfer??
lyapunov(K::SolutionInvWrapper,ac=acim(K)) = lyapunov(markovmap(Transfer(K)),ac)
lyapunov(L::Operator,ac=acim(L)) = lyapunov(markovmap(L),ac)

function lyapunov(m::MarkovMap,ac=acim(m)::Fun)
  lyap = zero(eltype(ac))
  for b in branches(m)
    lyap += sum(Fun(x->log(abs(mapD(b,convert(eltype(ac),x))))*ac(x),b.domain)) #TODO: why convert?
  end
  lyap
end
lyapunov(m::AbstractCircleMap,ac=acim(m)) =
  sum(Fun(x->log(abs(mapD(m,convert(eltype(ac),x))))*ac(x),m.domain))


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
(bi::MarkovBranchInverse)(x::Number) = mapinv(bi.b,x)
inv(b::MarkovBranch) = MarkovBranchInverse(b)
## TODO one day: convert routine to MarkovBranch?

immutable MarkovBranchDerivativeInverse{B<:MarkovBranch}
  b::B
end
(bi::MarkovBranchDerivativeInverse)(x::Number) = mapinvD(bd.b,x)
Base.ctranspose(bi::MarkovBranchInverse) = MarkovBranchDerivativeInverse(bi.b)
# Base.transpose(b::MarkovBranch) = MarkovBranchDerivative(b)

immutable MarkovMapDerivative{M<:AbstractMarkovMap}
  m::M
end
(md::MarkovMapDerivative)(x::Number) = mapD(md.m,x)
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
