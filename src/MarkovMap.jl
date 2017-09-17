# MarkovMaps
export MarkovMap, branch, nbranches, modulomap, induce, CircleMap

@compat abstract type AbstractMarkovMap{D<:Domain,R<:Domain} end# <: Function
#abstract AbstractDerivativeMarkovMap{D<:Domain,R<:Domain,T,FF} <: AbstractMarkovMap{D,R,T,FF}

Base.summary(m::AbstractMarkovMap) =  string(typeof(m).name.name)*":"*string(domain(m))*"↦"*string(rangedomain(m)) #branches??
# Base.eltype(m::AbstractMarkovMap) = eltype(rangedomain(m))
# Base.show(io::IO,m::AbstractMarkovMap) = print(io,typeof(m)) #temporary


@compat struct MarkovMap{D<:Domain,R<:Domain,B<:MarkovBranch} <: AbstractMarkovMap{D,R}
  branches::AbstractVector{B}
  domain::D
  rangedomain::R
  @compat function MarkovMap{D,R,B}(branches::AbstractVector{B},dom::D,ran::R) where
          {B<:MarkovBranch, D<:Domain, R<:Domain}
    @assert all(b->issubset(b.domain,dom),branches)
    @assert all(b->(b.rangedomain == ran),branches)
    # for i = 1:length(branches)
    #   for j = 1:i-1
    #     ~isempty(branches[i].domain ∩ branches[j].domain) && error("Overlapping domains in branches $i and $j")
    #   end
    # end
    # length(dom) == 1 && arclength(dom)/sum([arclength(b.domain) for b in branches]) < 1-200eps(eltype(dom)) && warn("Warning: possibly missing branches")
    new(branches,dom,ran)
  end
end
MarkovMap(branches::AbstractVector,dom,ran) = MarkovMap(branches,Domain(dom),Domain(ran))

function MarkovMap{B<:MarkovBranch}(branches::AbstractVector{B},dom,ran)
  domd = Domain(dom); randm = Domain(ran)
  MarkovMap{typeof(domd),typeof(randm),B}(branches,domd,randm)
end

coveringsegment(ds::AbstractArray) = coveringsegment([Domain(d) for d in ds])
coveringsegment{T<:Domain}(dsm::AbstractArray{T}) = Segment(minimum(first(Domain(d)) for d in dsm),maximum(last(Domain(d)) for d in dsm))

function MarkovMap(fs::AbstractVector,ds::AbstractVector,ran=coveringsegment(ds);dir=Forward,
          diff=[autodiff(fs[i],(dir==Forward ? ds[i] : ran)) for i in eachindex(fs)])
  @assert length(fs) == length(ds)
  randm=Domain(ran)
  MarkovMap([branch(fs[i],Domain(ds[i]),randm,diff[i];dir=dir,ftype=eltype(fs),difftype=eltype(diff)) for i in eachindex(fs)],
        coveringsegment(ds),Domain(ran))
end

(m::MarkovMap)(i::Integer,x) = (m.branches[i])(x)
(m::MarkovMap)(x) = (m.branches[getbranch(m,x)])(x)
for FUN in (:mapD,:mapP)
 @eval $FUN(m::MarkovMap,x) = $FUN(m.branches[getbranch(m,x)],x)
end
for FUN in (:mapD,:mapP,:mapinv,:mapinvD,:mapinvP)
  @eval $FUN(m::MarkovMap,i::Integer,x) = $FUN(m.branches[i],x)
end

branches(m) = m.branches
nbranches(m::MarkovMap) = length(m.branches)
nneutral(m::MarkovMap) = sum([isa(b,NeutralBranch) for b in m.branches])
getbranch(m::MarkovMap,x) = in(x,m.domain) ? findfirst([in(x,domain(b)) for b in m.branches]) : throw(DomainError)

# # TODO: maybe roll into branch constructors? maybe remove??
# function MarkovMap{D,R,ff}(
#     f::AbstractVector{ff},dfdx::AbstractVector{ff},dom::D,ran::R)
#   T = eltype(R)
#   bl = Array{T}(length(f)); bu = Array{T}(length(f))
#   tol = 10eps(maximum(abs([dom.a,dom.b])))
#   for i = 1:length(f)
#     ba = interval_newton(f[i],dfdx[i],ran.a,dom.a,dom.b,interval_guess(x,dom,ran),tol)
#     bb = interval_newton(f[i],dfdx[i],ran.b,dom.a,dom.b,interval_guess(x,dom,ran),tol)
#     bl[i] = min(ba,bb); bu[i] = max(ba,bb)
#     bl[i] < dom.a && bl[i]-dom.a > -2.5tol && (bl[i] = dom.a)
#     bu[i] > dom.b && bu[i]-dom.b < 2.5tol && (bu[i] = dom.b)
#   end
#   for i = 1:length(f)
#     for j = 1:length(f)
#       abs(bu[j]-bl[i]) < 2.5tol && ((j ==i) ?
#                                       error("Width of branch $i too small for automated branch edge checking...") :
#                                       (bu[j] = bl[i] = (bu[j]+bl[i])/2))
#     end
#   end
#   MarkovMap(f,dfdx,dom,ran,bl,bu)
# end
# #MarkovMap(dom::D,ran::R,f::AbstractVector{ff},dfdx::AbstractVector{ff}) = MarkovMap{D,R,eltype(R),ff}(dom,ran,f,dfdx)
# MarkovMap{ff<:ApproxFun.Fun}(f::AbstractVector{ff},dom::Domain,ran::Domain) = MarkovMap(f,[fi' for fi in f],dom,ran)
#
# # function MarkovMap{ff,gg}(dom::Domain,ran::Domain,f::AbstractVector{ff},dfdx::AbstractVector{gg},args...;kwargs...)
# #   pp = promote_type(ff,gg)
# #   MarkovMap(dom,ran,convert(Vector{pp},f),convert(Vector{pp},dfdx),args...;kwargs...)
# # end
# MarkovMap{ff}(f::AbstractVector{ff},dom::Domain,args...;kwargs...) = MarkovMap(f,dom,dom,args...;kwargs...)

# nice constructors

@compat struct Offset{F,T}
  f::F
  offset::T
end
(of::Offset)(x) = of.f(x)-of.offset

function modulomap{ff}(f::ff,dom,ran=dom;diff=autodiff(f,dom))
  domd = Domain(dom); randm = Domain(ran)
  fa = f(first(domd)); fb = f(last(domd))
  L = arclength(randm)
  nb_est = (fb-fa)/L; nb = round(Int,nb_est) # number of branches
  σ = sign(nb); NB = abs(nb)

  # @assert in(fa,∂(randm))
  # @assert nb_est ≈ nb
  @assert nb != 0

  breakpoints = Array{eltype(domd)}(NB+1)
  breakpoints[1] = first(domd)

  for i = 1:NB-1
    breakpoints[i+1] = domain_newton(f,diff,fa+σ*i*L,domd)
  end
  breakpoints[end] = last(domd)

  fs = [Offset(f,(i-1)*σ*L) for i = 1:NB]
  ds = [Interval(breakpoints[i],breakpoints[i+1]) for i = 1:NB]
  MarkovMap(fs,ds,randm,diff=fill(diff,NB))
end







# Inducing

@compat struct InducedMarkovMap{M<:MarkovMap,B<:MarkovBranch,D<:Domain,R<:Domain} <: AbstractMarkovMap{D,R}
  m::M
  b::B
  domain::D
  rangedomain::R
  function InducedMarkovMap{M,B,D,R}(m::M,b::B,dom::D,ran::R) where {M<:MarkovMap,B<:MarkovBranch,D<:Domain,R<:Domain} #requires m,b: their domains -> m ∪ b
    @assert dom == ran
    @assert domain(m) == dom
    @assert rangedomain(b) == rangedomain(m)
    @assert domain(m) == setdiff(rangedomain(m),domain(b))
    new(m,b,dom,ran)
  end
end
InducedMarkovMap(m::MarkovMap,b::MarkovBranch,dom::Domain,ran::Domain) =
  InducedMarkovMap{typeof(m),typeof(b),typeof(dom),typeof(ran)}(m,b,dom,ran)

function mapinduceP(b::MarkovBranch,x) #hack
  y = copy(x)
  dy = one(y)
  while in(y,b.domain)
    p = unsafe_mapP(b,y)
    y == p[1] && error("Map cannot be induced at this point as a fixed point is in its forward orbit")
    y = p[1]
    dy *= p[2]
  end
  (y,dy)
end
mapinduce(b::MarkovBranch,x) = mapinduceP(b,x)[1]
mapinduceD(b::MarkovBranch,x) = mapinduceP(b,x)[2]


function induce(m::MarkovMap,n::Integer)
  @assert domain(m) == rangedomain(m)
  b = m.branches[n]
  dom = setdiff(m.domain,b.domain)
  InducedMarkovMap(MarkovMap(m.branches[setdiff(1:nbranches(m),n)],dom,m.rangedomain),b,dom,dom)
end
nneutral(m::InducedMarkovMap) = nneutral(m.m) # if your maps stop being uniformly expanding this is bad

function mapP(m::InducedMarkovMap,x)
  y = mapP(m.m,x)
  if in(y[1],domain(m.b))
    a = mapinduceP(m.b,y[1])
    y = [a[1],y[2]*a[2]]
  end
  y
end
(m::InducedMarkovMap)(x) = mapP(m,x)[1]
mapD(m::InducedMarkovMap,x) = mapP(m,x)[2]




# Domain calls

for TYP in (:MarkovMap,:InducedMarkovMap)
  @eval ApproxFun.domain(m::$TYP) = m.domain
  @eval rangedomain(m::$TYP) = m.rangedomain
end

# InverseCache

# @compat mutable struct MarkovInverseCache{M<:AbstractMarkovMap,S<:Space,T} <: AbstractMarkovMap
#   m::M
#   sp::S
#   cache::Dict{Int,Array{T,2}}
# end
# function MarkovInverseCache(m::AbstractMarkovMap,s::Space)
#   @assert rangedomain(m) == domain(s)
#   MarkovInverseCache{typeof(m),typeof(s),eltype(domain(m))}(m,s,Dict{Int,Array{eltype(domain(m)),2}}())
# end
#
# MarkovInverseCache(m::AbstractMarkovMap) = MarkovInverseCache(m,ApproxFun.Space(rangedomain(m)))
#
# ApproxFun.domain(mc::MarkovInverseCache) = domain(mc.m)
# rangedomain(mc::MarkovInverseCache) = rangedomain(mc.m)
# length(mc::MarkovInverseCache) = length(mc.m)
# rangespace(mc::MarkovInverseCache) = mc.sp
# function extenddata!(m::MarkovInverseCache,n::Int)
#   v = Array{eltype(m)}(length(m),n)
#   pts = points(rangespace(m),n)
#   for k = 1:n
#     for i = 1:length(m)
#       v[i,k] = mapinv(m.m,i,pts[k])
#     end
#   end
#   push!(m.cache,n=>v)
# end
#
# for FN in (:mapinv,:mapinvD)
#   @eval $FN(m::MarkovInverseCache,i::Integer,x) = $FN(m.m,i,x)
# end
#
# function mapinv(m::MarkovInverseCache,i::Integer,x::InterpolationNode)
#   rangespace(m) != space(x) && return mapinv(m,i,convert(eltype(x),x))
#   ~haskey(m.cache,x.n) && extenddata!(m,x.n)
#   m.cache[x.n][i,x.k]
# end
# # mapinvD{M<:InverseDerivativeMarkovMap}(m::MarkovInverseCache{M},i::Integer,x::InterpolationNode) =
# #   mapinvD(m.m,i,x)
# mapinvD(m::MarkovInverseCache,i::Integer,x::InterpolationNode) =
#   isa(m.branches[i],RevExpandingBranch) ? mapinvD(m.m,i,mapinv(m,i,x)) : inv(mapD(m.m,i,mapinv(m,i,x)))
# # mapDsign(m::MarkovInverseCache,i::Integer) = mapDsign(m.m,i)
