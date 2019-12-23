# MarkovMaps
export MarkovMap, IntervalMap, branch, nbranches, modulomap, CircleMap

@compat abstract type AbstractIntervalMap{D<:Domain,R<:Domain} end# <: Function
@compat abstract type AbstractMarkovMap{D<:Domain,R<:Domain} <: AbstractIntervalMap{D,R} end# <: Function
# @compat abstract type AbstractMarkovMap{D<:Domain,R<:Domain} end# <: Function
#abstract AbstractDerivativeMarkovMap{D<:Domain,R<:Domain,T,FF} <: AbstractMarkovMap{D,R,T,FF}

Base.show(io::IO, m::AbstractIntervalMap) = print(io, string(typeof(m).name.name)*
    " "*string(domain(m))*" → "*string(rangedomain(m))*" with $(nbranches(m)) branches") #branches??
prectype(m::AbstractIntervalMap) = promote_type(prectype(domain(m)),prectype(rangedomain(m)))
eltype(m::AbstractMarkovMap) = promote_type(eltype(domain(m)),eltype(rangedomain(m)))

# Base.show(io::IO,m::AbstractMarkovMap) = print(io,typeof(m)) #temporary

# Domain calls

ApproxFun.domain(m::AbstractIntervalMap) = m.domain
rangedomain(m::AbstractIntervalMap) = m.rangedomain

# suppose you can only put nfps in a map, in general (<== What does this mean??)
containsnfp(d,m) = any(neutralfixedpoints(m) .∈ d)
neutralfixedpoints(m::AbstractIntervalMap) = 0


"""
    MarkovMap(branches::Vector, domain, rangedomain)

Generate a computer representation of a full-branch uniformly-expanding interval map `domain` → `rangedomain` using a vector describing the branches of the map.
"""
@compat struct MarkovMap{D<:Domain,R<:Domain,B<:ExpandingBranch} <: AbstractMarkovMap{D,R}
  branches::Vector{B}
  domain::D
  rangedomain::R
  @compat function MarkovMap{D,R,B}(branches::AbstractVector{B},dom::D,ran::R) where
          {B<:ExpandingBranch, D<:Domain, R<:Domain}
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
MarkovMap(branches::Vector{B},dom::D,ran::R) where {D<:Domain,R<:Domain,B<:ExpandingBranch}= MarkovMap{D,R,B}(branches,dom,ran)
MarkovMap(branches::Vector{B},dom,ran) where {B<:ExpandingBranch} = MarkovMap(branches,convert(Domain,dom),convert(Domain,ran))

function MarkovMap(branches::AbstractVector{B},dom,ran) where {B<:ExpandingBranch}
  domd = convert(Domain,dom); randm = convert(Domain,ran)
  MarkovMap{typeof(domd),typeof(randm),B}(branches,domd,randm)
end

"""
    MarkovMap(fs::Vector, ds::Vector, ran = coveringinterval(ds); dir=Forward, diff=autodiff(...)))

Generate a MarkovMap with branches given by elements of `fs`` defined on subdomains given by `ds`, onto a vector `ran`.

The keyword argument `dir` stipulates whether the elements of `fs` are the branches (`Forward`) or the branches' inverses (`Reverse`).

The keyword argument `diff` provides the derivatives of the `fs`. By default it is the automatic derivatives of `fs`.
"""
function MarkovMap(fs::AbstractVector,ds::AbstractVector,ran=coveringinterval(ds);dir=Forward,
          diff=[autodiff(fs[i],(dir==Forward ? ds[i] : ran)) for i in eachindex(fs)])
  @assert length(fs) == length(ds)
  randm=convert(Domain,ran)
  MarkovMap([branch(fs[i],convert(Domain,ds[i]),randm,diff[i];dir=dir,ftype=eltype(fs),difftype=eltype(diff)) for i in eachindex(fs)],
        coveringinterval(ds),convert(Domain,ran))
end


# TODO: MarkovMap not requiring branch entries?

struct IntervalMap{D<:Domain,R<:Domain,B<:AbstractBranch} <: AbstractIntervalMap{D,R}
  branches::Vector{B}
  domain::D
  rangedomain::R
  @compat function IntervalMap{D,R,B}(branches::Vector{B},dom::D,ran::R) where
          {B<:ExpandingBranch, D<:Domain, R<:Domain}
    @assert all(b->approx_issubset(b.domain,dom),branches)
    @assert all(b->approx_issubset(b.rangedomain,ran),branches)
    new(branches,dom,ran)
  end
end
IntervalMap(branches::Vector{B},dom::Domain,ran::Domain) where {B<:AbstractBranch} = IntervalMap{typeof(dom),typeof(ran),eltype(branches)}(branches,dom,ran)
IntervalMap(branches::Vector{B},dom,ran) where {B<:AbstractBranch} = IntervalMap(branches,convert(Domain,dom),convert(Domain,ran))

function IntervalMap(branches::AbstractVector{B},dom,ran) where {B<:ExpandingBranch}
  domd = convert(Domain,dom); randm = convert(Domain,ran)
  IntervalMap{typeof(domd),typeof(randm),B}(branches,domd,randm)
end

function IntervalMap(fs::AbstractVector,ds::AbstractVector,
          ran=coveringinterval([mapinterval(fs[i],ds[i]) for i in eachindex(fs)]);dir=Forward,
          diff=[autodiff(fs[i],(dir==Forward ? ds[i] : ran)) for i in eachindex(fs)])
  rand = convert(Domain,ran)
  @assert length(fs) == length(ds)
  # println(rand)
  # println(extrema(mapinterval(fs[1],ds[1])))
  # println((approx_issubset(mapinterval(fs[i],ds[i]),rand) for i in eachindex(fs))|>collect)
  @assert all(approx_issubset(mapinterval(fs[i],ds[i]),rand) for i in eachindex(fs))
  IntervalMap([branch(fs[i],convert(Domain,ds[i]),mapinterval(fs[i],ds[i]),diff[i];dir=dir,ftype=eltype(fs),difftype=eltype(diff)) for i in eachindex(fs)],
        coveringinterval(ds),rand)
end

const SimpleBranchedMap{D<:Domain,R<:Domain,B<:AbstractBranch} = Union{MarkovMap{D,R,B},IntervalMap{D,R,B}} #?? ComposedMap

for FUN in (:mapD,:mapP)
 @eval $FUN(m::SimpleBranchedMap,x) = $FUN(m.branches[getbranchind(m,x)],x)
end
for FUN in (:mapD,:mapP,:mapinv,:mapinvD,:mapinvP)
  @eval $FUN(m::SimpleBranchedMap,i::Integer,x) = $FUN(m.branches[i],x)
end
for TYP = (:MarkovMap, :IntervalMap)
  @eval begin
    (m::$TYP)(i::Integer,x) = (m.branches[i])(x)
    (m::$TYP)(x) = (m.branches[getbranchind(m,x)])(x)
  end
end

"""
    branches(m)

Return the branches of the map `m`.
"""
branches(m::SimpleBranchedMap) = m.branches

"""
    nbranches(m)

Return the number of branches of `m`.
"""
nbranches(m::SimpleBranchedMap) = length(m.branches)

"""
    eachbranchindex(m)

Return an iterator giving the indices of the branches of `m`
"""
eachbranchindex(m::SimpleBranchedMap) = 1:nbranches(m)

neutralfixedpoints(m::SimpleBranchedMap) = [nfp(b) for b in branches(m)[isa.(NeutralBranch,branches(m))]]
nneutral(m::SimpleBranchedMap) = sum([isa(b,NeutralBranch) for b in m.branches])
getbranchind(m::SimpleBranchedMap,x) = x ∈ domain(m) ? findfirst([(x ∈ domain(b)) for b in m.branches]) : error("DomainError: $x ∉ $(domain(m))")
getbranchind(m::SimpleBranchedMap,x::AbstractInterval) = (x ⊆ domain(m)) ? findfirst([(x ⊆ domain(b)) for b in m.branches]) : error("DomainError: $x ⊈ $(domain(m))")
getbranchind(m::SimpleBranchedMap,x::Interval) = getbranchind(m,convert(Domain,x))
branchindtype(m::SimpleBranchedMap) = Int

# Transfer function

function transferfunction(x,m::SimpleBranchedMap,f)
  y = zero(eltype(x));
  for b in branches(m)
    y += transferbranch(x,b,f)
  end;
  y
end

function transferfunction_int(x,y,m::SimpleBranchedMap,sk)
  q = zero(eltype(x));

  for b in branches(m)
    q += transferbranch_int(x,y,b,sk)
  end;
  q
end

# nice constructors

# modulomap

@compat struct FwdOffset{F,T}
  f::F
  offset::T
end
(of::FwdOffset)(x) = of.f(x)-of.offset

function forwardmodulomap(f::ff,dom,ran=dom,diff=autodiff(f,dom)) where {ff}
  domd = convert(Domain,dom); randm = convert(Domain,ran)
  fa = f(leftendpoint(domd)); fb = f(rightendpoint(domd))
  L = arclength(randm)
  nb_est = (fb-fa)/L; nb = round(Int,nb_est) # number of branches
  σ = sign(nb); NB = abs(nb)

  # @assert in(fa,∂(randm))
  # @assert nb_est ≈ nb
  @assert nb != 0

  breakpoints = Array{eltype(domd)}(undef,NB+1)
  breakpoints[1] = leftendpoint(domd)

  for i = 1:NB-1
    breakpoints[i+1] = domain_newton(f,diff,fa+σ*i*L,domd)
  end
  breakpoints[end] = rightendpoint(domd)

  fs = [FwdOffset(f,(i-(σ == 1))*σ*L) for i = 1:NB]
  ds = [Interval(breakpoints[i],breakpoints[i+1]) for i = 1:NB]
  MarkovMap(fs,ds,randm,diff=fill(diff,NB))
end


@compat struct RevOffset{F,T}
  fn::F
  offset::T
end
(ov::RevOffset)(x) = ov.fn(x+ov.offset)

reversemodulomap(args...) = _NI("reversemodulomap")
# TODO
# function reversemodulomap{ff}(v::ff,dom,ran=dom,diff=autodiff(v,dom);maxcover=10000)
#   domd = convert(Domain,dom); randm = convert(Domain,ran)
#   dr = arclength(randm); dd = arclength(domd)
#   ra = leftendpoint(randm); va = v(ra); vb = v()
#   in(va,∂(domd)) || error("v($a) not in boundary of domain")
#   flip = (va ≈ rightendpoint(domd))
#   drflip = (-1)^flip * dr
#
#   NB = 1; vr = v(ra+drflip)-va
#   while (abs(vr) <= dd && NB < maxcover)
#     vr = v(ra+NB*drflip)-va
#     abs(vr) ≈ dd && break
#     abs(vr) > dd && error("Inverse lift doesn't appear to have an inverse")
#     NB += 1
#   end
#   NB == maxcover && error("Can't get to the end of the inverse lift after 10000 steps")
#
#   rng = flip ? 0:NB : NB:-1:1
#   breakpoints = v.(ra+rng*dr)
#
#   vs = [RevOffset(v,ra+(i-1)*dr) for i = rng]
#   dvdxs = [RevOffset(diff,ra+(i-1)*dr) for i = rng]
#   ds = [Interval(breakpoints[i],breakpoints[i+1]) for i = 1:NB]
#   MarkovMap(vs,ds,randm,diff=dvdxs,dir=Reverse)
# end

"""
    modulomap(f, D, R=dom; diff= autodiff(f,dom))

Output MarkovMap or CircleMap m: D → R such that m(x) = f(x) mod R.
"""
modulomap(f,d,r=d;dir=Forward,diff=autodiff(f,dir==Forward ? d : r)) =
    dir == Forward ? forwardmodulomap(f,d,r,diff) : reversemodulomap(f,d,r,diff)
