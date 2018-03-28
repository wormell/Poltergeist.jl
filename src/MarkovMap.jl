# MarkovMaps
export MarkovMap, IntervalMap, branch, nbranches, modulomap, induce, CircleMap

@compat abstract type AbstractIntervalMap{D<:Domain,R<:Domain} end# <: Function
@compat abstract type AbstractMarkovMap{D<:Domain,R<:Domain} <: AbstractIntervalMap{D,R} end# <: Function
# @compat abstract type AbstractMarkovMap{D<:Domain,R<:Domain} end# <: Function
#abstract AbstractDerivativeMarkovMap{D<:Domain,R<:Domain,T,FF} <: AbstractMarkovMap{D,R,T,FF}

Base.summary(m::AbstractIntervalMap) =  string(typeof(m).name.name)*":"*string(domain(m))*"↦"*string(rangedomain(m))*" with $(nbranches(m)) branches" #branches??
# Base.eltype(m::AbstractMarkovMap) = eltype(rangedomain(m))
# Base.show(io::IO,m::AbstractMarkovMap) = print(io,typeof(m)) #temporary


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
MarkovMap{D<:Domain,R<:Domain,B<:ExpandingBranch}(branches::Vector{B},dom::D,ran::R) = MarkovMap{D,R,B}(branches,dom,ran)
MarkovMap{B<:ExpandingBranch}(branches::Vector{B},dom,ran) = MarkovMap(branches,Domain(dom),Domain(ran))

function MarkovMap{B<:ExpandingBranch}(branches::AbstractVector{B},dom,ran)
  domd = Domain(dom); randm = Domain(ran)
  MarkovMap{typeof(domd),typeof(randm),B}(branches,domd,randm)
end

function MarkovMap(fs::AbstractVector,ds::AbstractVector,ran=mapinterval(first(fs),first(ds));dir=Forward,
          diff=[autodiff(fs[i],(dir==Forward ? ds[i] : ran)) for i in eachindex(fs)])
  @assert length(fs) == length(ds)
  randm=Domain(ran)
  MarkovMap([branch(fs[i],Domain(ds[i]),randm,diff[i];dir=dir,ftype=eltype(fs),difftype=eltype(diff)) for i in eachindex(fs)],
        coveringsegment(ds),Domain(ran))
end


# TODO: MarkovMap not requiring branch entries?

type IntervalMap{D<:Domain,R<:Domain,B<:AbstractBranch} <: AbstractIntervalMap{D,R}
  branches::Vector{B}
  domain::D
  rangedomain::R
  @compat function IntervalMap{D,R,B}(branches::Vector{B},dom::D,ran::R) where
          {B<:ExpandingBranch, D<:Domain, R<:Domain}
    @assert all(b->issubset(b.domain,dom),branches)
    @assert all(b->issubset(b.rangedomain,ran),branches)
    new(branches,dom,ran)
  end
end
IntervalMap{B<:AbstractBranch}(branches::Vector{B},dom::Domain,ran::Domain) = IntervalMap{typeof(dom),typeof(ran),eltype(branches)}(branches,dom,ran)
IntervalMap{B<:AbstractBranch}(branches::Vector{B},dom,ran) = IntervalMap(branches,Domain(dom),Domain(ran))

function IntervalMap{B<:ExpandingBranch}(branches::AbstractVector{B},dom,ran)
  domd = Domain(dom); randm = Domain(ran)
  IntervalMap{typeof(domd),typeof(randm),B}(branches,domd,randm)
end

function IntervalMap(fs::AbstractVector,ds::AbstractVector,
          ran=coveringsegment([mapinterval(fs[i],ds[i]) for i in eachindex(fs)]);dir=Forward,
          diff=[autodiff(fs[i],(dir==Forward ? ds[i] : ran)) for i in eachindex(fs)])
  @assert length(fs) == length(ds)
  randm=Domain(ran)
  IntervalMap([branch(fs[i],Domain(ds[i]),randm,diff[i];dir=dir,ftype=eltype(fs),difftype=eltype(diff)) for i in eachindex(fs)],
        coveringsegment(ds),Domain(ran))
end

const BranchedMap{D<:Domain,R<:Domain} = Union{MarkovMap{D,R},IntervalMap{D,R}} #?? ComposedMap

for TYP = (:MarkovMap, :IntervalMap)
  for FUN in (:mapD,:mapP)
   @eval $FUN(m::$TYP,x) = $FUN(m.branches[getbranch(m,x)],x)
  end
  for FUN in (:mapD,:mapP,:mapinv,:mapinvD,:mapinvP)
    @eval $FUN(m::$TYP,i::Integer,x) = $FUN(m.branches[i],x)
  end
  @eval begin
    (m::$TYP)(i::Integer,x) = (m.branches[i])(x)
    (m::$TYP)(x) = (m.branches[getbranch(m,x)])(x)

    branches(m::$TYP) = m.branches
    nbranches(m::$TYP) = length(m.branches)
    eachbranchindex(m::$TYP) = 1:nbranches(m)

    nneutral(m::$TYP) = sum([isa(b,NeutralBranch) for b in m.branches])
    getbranch(m::$TYP,x) = temp_in(x,m.domain) ? findfirst([temp_in(x,domain(b)) for b in m.branches]) : error("DomainError: $x ∉ $(m.domain)")

    # Domain calls

    ApproxFun.domain(m::$TYP) = m.domain
    rangedomain(m::$TYP) = m.rangedomain

    # Transfer function

    function transferfunction(x,m::$TYP,f,T)
      y = zero(eltype(x));
      for b in branches(m)
        y += transferbranch(x,b,f,T)
      end;
      y
    end

    function transferfunction_int(x,y,m::$TYP,sk,T)
      q = zero(eltype(x));

      for b in branches(m)
        q += transferbranch_int(x,y,b,sk,T)
      end;
      q
    end
  end
end

# nice constructors

# modulomap

  @compat struct FwdOffset{F,T}
    f::F
    offset::T
  end
  (of::FwdOffset)(x) = of.f(x)-of.offset

function forwardmodulomap{ff}(f::ff,dom,ran=dom,diff=autodiff(f,dom))
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

  fs = [FwdOffset(f,(i-(σ == 1))*σ*L) for i = 1:NB]
  ds = [Interval(breakpoints[i],breakpoints[i+1]) for i = 1:NB]
  MarkovMap(fs,ds,randm,diff=fill(diff,NB))
end


@compat struct RevOffset{F,T}
  fn::F
  offset::T
end
(ov::RevOffset)(x) = ov.fn(x+ov.offset)

# TODO
# function reversemodulomap{ff}(v::ff,dom,ran=dom,diff=autodiff(v,dom);maxcover=10000)
#   domd = Domain(dom); randm = Domain(ran)
#   dr = arclength(randm); dd = arclength(domd)
#   ra = first(randm); va = v(ra); vb = v()
#   in(va,∂(domd)) || error("v($a) not in boundary of domain")
#   flip = (va ≈ last(domd))
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

modulomap(f,d,r=d;dir=Forward,diff=autodiff(f,dir==Forward ? d : r)) =
    dir == Forward ? forwardmodulomap(f,d,r,diff) : reversemodulomap(f,d,r,diff)
