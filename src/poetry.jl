export perturb

# Markov derivatives
@compat struct ExpandingBranchDerivative{B<:ExpandingBranch}
  b::B
end
(bd::ExpandingBranchDerivative)(x::Number) = mapD(bd.b,x)
Base.ctranspose(b::ExpandingBranch) = ExpandingBranchDerivative(b)
# Base.transpose(b::ExpandingBranch) = ExpandingBranchDerivative(b)

@compat struct ExpandingBranchInverse{B<:ExpandingBranch}
  b::B
end
(bi::ExpandingBranchInverse)(x::Number) = mapinv(bi.b,x)
inv(b::ExpandingBranch) = ExpandingBranchInverse(b)
## TODO one day: convert routine to ExpandingBranch?

@compat struct ExpandingBranchDerivativeInverse{B<:ExpandingBranch}
  b::B
end
(bi::ExpandingBranchDerivativeInverse)(x::Number) = mapinvD(bd.b,x)
Base.ctranspose(bi::ExpandingBranchInverse) = ExpandingBranchDerivativeInverse(bi.b)
# Base.transpose(b::ExpandingBranch) = ExpandingBranchDerivative(b)

@compat struct MarkovMapDerivative{M<:AbstractMarkovMap}
  m::M
end
(md::MarkovMapDerivative)(x::Number) = mapD(md.m,x)
Base.ctranspose(b::AbstractMarkovMap) = MarkovMapDerivative(b)

# Linear response perturbations
perturb(d,X,ϵ) = perturb(Domain(d),X,ϵ)
perturb(d::IntervalDomain,X,ϵ) = MarkovMap([x->x+ϵ*X(x)],[d],d)
perturb(d::PeriodicDomain,X,ϵ) = FwdCircleMap([x->x+ϵ*X(x)],d)
perturb(m::AbstractMarkovMap,X,ϵ) = perturb(rangedomain(m),X,ϵ)∘m

# Eigvals overloads
eigvals(m::AbstractMarkovMap,n) = eigvals(Transfer(m),n)

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
