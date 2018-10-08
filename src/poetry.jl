export perturb

# Markov derivatives
@compat struct ExpandingBranchDerivative{B<:ExpandingBranch}
  b::B
end
(bd::ExpandingBranchDerivative)(x::Number) = mapD(bd.b,x)
Base.adjoint(b::ExpandingBranch) = ExpandingBranchDerivative(b)
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
Base.adjoint(bi::ExpandingBranchInverse) = ExpandingBranchDerivativeInverse(bi.b)
# Base.transpose(b::ExpandingBranch) = ExpandingBranchDerivative(b)

@compat struct MarkovMapDerivative{M<:AbstractMarkovMap}
  m::M
end
(md::MarkovMapDerivative)(x::Number) = mapD(md.m,x)
Base.adjoint(b::AbstractMarkovMap) = MarkovMapDerivative(b)

# Linear response perturbations
"""
    perturb(d, X, ϵ)

Construct a self-map on domain d: x ↦ x + ϵ X(x)
"""
perturb(d,X,ϵ) = perturb(convert(Domain,d),X,ϵ)
perturb(d::IntervalDomain,X,ϵ) = MarkovMap([x->x+ϵ*X(x)],[d],d)
perturb(d::PeriodicDomain,X,ϵ) = FwdCircleMap([x->x+ϵ*X(x)],d)

"""
    perturb(m::AbstractMarkovMap, X, ϵ)

Output perturbation of m: x ↦ m(x) + ϵ X(m(x))
"""
perturb(m::AbstractMarkovMap,X,ϵ) = perturb(rangedomain(m),X,ϵ)∘m

# Eigvals overloads
"""
    eigvals(m::AbstractMarkovMap, n)

Output eigenvalues of Transfer(m) using n×n Galerkin discretisation.

Calls directly to ApproxFun: you can also call eigvals(Transfer(m), n)
"""
@compat LinearAlgebra.eigvals(m::AbstractMarkovMap,n::Int64) = LinearAlgebra.eigvals(Transfer(m),n)
ApproxFun.eigs(m::AbstractMarkovMap,n::Int64) = ApproxFun.eigs(Transfer(m),n)

"""
    eigvecs(m::AbstractMarkovMap, n)

Output eigenfunctions of Transfer(m) using n×n Galerkin discretisation.

Calls directly to ApproxFun: you can also call eigvecs(Transfer(m), n)
"""
@compat LinearAlgebra.eigvecs(m::AbstractMarkovMap,n::Int64) = ApproxFun.eigvecs(Transfer(m),n)

# # plotting
# function plot(m::MarkovMap)
#   pts = cfstype(m)[]
#   vals = cfstype(m)[]
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
