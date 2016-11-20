export MarkovMap, branch, nbranches

# MarkovBranch

abstract MarkovBranch{D<:Domain,R<:Domain}

Base.summary(b::MarkovBranch) =  string(typeof(b).name.name)*":"*string(domain(b))*"↦"*string(rangedomain(b)) #branches??
Base.eltype(b::MarkovBranch) = eltype(rangedomain(b))
Base.show(io::IO,b::MarkovBranch) = print(io,typeof(b)) #temporary

immutable NeutralBranch{D,R} <: MarkovBranch{D,R}
end

immutable FwdExpandingBranch{ff,gg,D<:Domain,R<:Domain} <: MarkovBranch{D,R}
  f::ff
  dfdx::gg
  domain::D
  rangedomain::R
  #  sgn::T
  function FwdExpandingBranch(fc,dfdxc,dom,ran)
    ~isempty(∂(ran)) && @assert all([in(fc(p),∂(ran)) for p in ∂(dom)])
    new(fc,dfdxc,dom,ran)
  end
end
FwdExpandingBranch{D<:Domain,R<:Domain,ff,gg}(f::ff,dfdx::gg,dom::D,ran::R) =
  FwdExpandingBranch{typeof(f),typeof(dfdx),typeof(dom),typeof(ran)}(f,dfdx,dom,ran)

unsafe_call(b::FwdExpandingBranch,x::Number) = b.f(x)
unsafe_mapD(b::FwdExpandingBranch,x::Number) = b.dfdx(x)
mapinv(b::FwdExpandingBranch,x::Number) = interval_newton(b.f,b.dfdx,x,b.domain)
mapinvD(b::FwdExpandingBranch,x::Number) = 1/b.dfdx(mapinv(b,x))

immutable RevExpandingBranch{ff,gg,D<:Domain,R<:Domain} <: MarkovBranch{D,R}
  v::ff
  dvdx::gg
  domain::D
  rangedomain::R
  #   sgn::T
  function RevExpandingBranch(vc,dvdxc,dom,ran)
    new(vc,dvdxc,dom,ran)
  end
end
RevExpandingBranch{ff,gg,D<:Domain,R<:Domain}(v::ff,dvdx::gg,dom::D,ran::R) =
  RevExpandingBranch{typeof(v),typeof(dvdx),typeof(dom),typeof(ran)}(v,dvdx,dom,ran)

unsafe_call(b::RevExpandingBranch,x::Number) = interval_newton(b.v,b.dvdx,x,b.rangedomain)
unsafe_mapD(b::RevExpandingBranch,x::Number) =  1/b.dvdx(unsafe_call(b,x))
mapinv(b::RevExpandingBranch,x::Number) = b.v(x)
mapinvD(b::RevExpandingBranch,x::Number) = b.dvdx(x)

for TYP in (:FwdExpandingBranch,:RevExpandingBranch)
  @eval (b::$TYP)(x::Number) = in(x,b.domain) ? unsafe_call(b,x) : throw(DomainError)
  @eval mapD(b::$TYP,x::Number) = in(x,b.domain) ? unsafe_mapD(b,x) : throw(DomainError)
  @eval domain(b::$TYP) = b.domain
  @eval rangedomain(B::$TYP) = b.rangedomain
end

# branch constructors

branch{ff,gg,D<:Domain,R<:Domain}(f::ff,dfdx::gg,dom::D,ran::R,dir::String="fwd",ftype::Type{ff}=ff,dfdxtype::Type{gg}=gg) =
  dir=="fwd" ? FwdExpandingBranch{ftype,dfdxtype,typeof(dom),typeof(ran)}(f,dfdx,dom,ran) : RevExpandingBranch{ftype,dfdxtype,typeof(dom),typeof(ran)}(f,dfdx,dom,ran)


function branch{ff,gg,D<:Domain,R<:Domain}(fs::Vector{ff},dfdxs::Vector{gg},doms::Vector{D},ran::R,dir::String="fwd")
  @assert length(fs) == length(dfdxs)
  @assert length(doms) == length(fs)
  ftype = promote_type([typeof(f) for f in fs]...)
  dfdxtype = promote_type([typeof(f) for f in dfdxs]...)
  TYP = dir=="fwd" ? FwdExpandingBranch{ftype,dfdxtype,D,R} : #eltype(doms),typeof(ran)
  RevExpandingBranch{ftype,dfdxtype,D,R}
  TYP[branch(fs[i],dfdxs[i],doms[i],ran,dir,ftype,dfdxtype) for i = 1:length(fs)]
end
function branch{ff,gg,T<:Number,R<:Domain}(fs::Vector{ff},dfdxs::Vector{gg},bl::Vector{T},bu::Vector{T},ran::R,dir::String="fwd")
  @assert length(bl) == length(fs)
  @assert length(bu) == length(fs)
  branch(fs,dfdxs,ApproxFun.Interval{T}[Interval(bl[i],bu[i]) for i = 1:length(bl)],ran,dir)
end

function branch{ff,gg,T<:Number,R<:Domain}(fs::Vector{ff},dfdxs::Vector{gg},b::Vector{T},ran::R,dir::String="fwd")
  @assert length(fs) == (length(b)-1)
  @assert issorted(b)
  branch(fs,dfdxs,ApproxFun.Interval{T}[Interval(b[i],b[i+1]) for i = 1:length(b)-1],ran,dir)
end

branch{ff,T<:Union{Number,Domain}}(fs::Vector{ff},b::Vector{T},args...) =
  branch(fs,[forwarddiff(f) for f in fs],b,args...)
branch{T<:Union{Number,Domain}}(f,b::T,args...) = branch(f,forwarddiff(f),b,args...)

branch{ff<:Fun,T<:Union{Number,Domain}}(fs::Vector{ff},b::Vector{T},args...) =
  branch(fs,[f' for f in fs],b,args...)
branch{T<:Union{Number,Domain}}(f::Fun,b::T,args...) = branch(f,f',b,args...)


# MarkovMaps

abstract AbstractMarkovMap{D<:Domain,R<:Domain}# <: Function
#abstract AbstractDerivativeMarkovMap{D<:Domain,R<:Domain,T,FF} <: AbstractMarkovMap{D,R,T,FF}

Base.summary(m::AbstractMarkovMap) =  string(typeof(m).name.name)*":"*string(domain(m))*"↦"*string(rangedomain(m)) #branches??
Base.eltype(m::AbstractMarkovMap) = eltype(rangedomain(m))
Base.show(io::IO,m::AbstractMarkovMap) = print(io,typeof(m)) #temporary


immutable MarkovMap{D<:Domain,R<:Domain,B<:MarkovBranch} <: AbstractMarkovMap{D,R}
  domain::D
  rangedomain::R
  branches::Vector{B}
  function MarkovMap(dom,ran,branches)
    @assert all([issubset(b.domain,dom) for b in branches])
    @assert all([b.rangedomain == ran for b in branches])
    for i = 1:length(branches)
      for j = 1:i-1
        ~isempty(branches[i].domain ∩ branches[j].domain) && error("Overlapping domains in branches $i and $j")
      end
    end
    arclength(dom)/sum([arclength(b.domain) for b in branches]) < nextfloat(one(eltype(dom)),-200) && warn("Warning: possibly missing branches")
    new{D,R,B}(dom,ran,branches)
  end
end
MarkovMap{D<:Domain,R<:Domain,B<:MarkovBranch}(dom::D,ran::R,branches::Vector{B}) = MarkovMap{typeof(dom),typeof(ran),B}(dom,ran,branches)

# MarkovMap{T,ff<:Number}(dom::Domain,ran::Domain,f::Vector{ff},dfdx::Vector{ff},bl::Vector{T},bu::Vector{T}) =
#   MarkovMap{T,ff}(dom,ran,f,dfdx,bl,bu)

function MarkovMap{D<:Domain,R<:Domain}(
    dom::D,ran::R,v1::Vector,v2::Vector,dir::String="fwd")
  MarkovMap(dom,ran,branch(v1,v2,ran,dir))
end
function MarkovMap{D<:Domain,R<:Domain}(
    dom::D,ran::R,v1::Vector,v2::Vector,v3::Vector,dir::String="fwd")
  MarkovMap(dom,ran,branch(v1,v2,v3,ran,dir))
end
function MarkovMap{D<:Domain,R<:Domain}(
    dom::D,ran::R,v1::Vector,v2::Vector,v3::Vector,v4::Vector,dir::String="fwd")
  MarkovMap(dom,ran,branch(v1,v2,v3,v4,ran,dir))
end


# TODO: maybe roll into branch constructors? maybe remove??
function MarkovMap{D,R,ff}(
    dom::D,ran::R,f::Vector{ff},dfdx::Vector{ff})
  T = eltype(R)
  bl = Array(T,length(f)); bu = Array(T,length(f))
  tol = 10eps(maxabs([dom.a,dom.b]))
  for i = 1:length(f)
    ba = interval_newton(f[i],dfdx[i],ran.a,dom.a,dom.b,tol)
    bb = interval_newton(f[i],dfdx[i],ran.b,dom.a,dom.b,tol)
    bl[i] = min(ba,bb); bu[i] = max(ba,bb)
    bl[i] < dom.a && bl[i]-dom.a > -2.5tol && (bl[i] = dom.a)
    bu[i] > dom.b && bu[i]-dom.b < 2.5tol && (bu[i] = dom.b)
  end
  for i = 1:length(f)
    for j = 1:length(f)
      abs(bu[j]-bl[i]) < 2.5tol && ((j ==i) ?
                                      error("Width of branch $i too small for automated branch edge checking...") :
                                      (bu[j] = bl[i] = (bu[j]+bl[i])/2))
    end
  end
  MarkovMap(dom,ran,f,dfdx,bl,bu)
end
#MarkovMap(dom::D,ran::R,f::Vector{ff},dfdx::Vector{ff}) = MarkovMap{D,R,eltype(R),ff}(dom,ran,f,dfdx)
MarkovMap{ff<:ApproxFun.Fun}(dom,ran,f::Vector{ff}) = MarkovMap(dom,ran,f,[fi' for fi in f])

# function MarkovMap{ff,gg}(dom::Domain,ran::Domain,f::Vector{ff},dfdx::Vector{gg},args...)
#   pp = promote_type(ff,gg)
#   MarkovMap(dom,ran,convert(Vector{pp},f),convert(Vector{pp},dfdx),args...)
# end
MarkovMap{ff}(dom::Domain,f::Vector{ff},args...) = MarkovMap(dom,dom,f,args...)

branches(m) = m.branches
nbranches(m::MarkovMap) = length(m.branches)
nneutral(m::MarkovMap) = sum([isa(b,NeutralBranch) for b in m.branches])
getbranch(m::MarkovMap,x::Number) = findfirst([in(x,b.domain) for b in m.branches])

(m::MarkovMap)(i::Integer,x) = (m.branches[i])(x)
(m::MarkovMap)(x::Number) = (m.branches[getbranch(m,x)])(x)
mapD(m::MarkovMap,i::Integer,x) = mapD(m.branches[i],x)
mapD(m::MarkovMap,x::Number) = mapD(m.branches[getbranch(m,x)],x)
mapinv(m::MarkovMap,i::Integer,x::Number) = mapinv(m.branches[i],x)
mapinvD(m::MarkovMap,i::Integer,x::Number) = mapinvD(m.branches[i],x)
#mapDsign(m::MarkovMap,i::Integer) = m.sgns[i]

#getindex(m::InverseDerivativeMarkovMap,:dvdx) = m.dvdx

ApproxFun.domain(m::MarkovMap) = m.domain
rangedomain(m::MarkovMap) = m.rangedomain

type MarkovInverseCache{M<:AbstractMarkovMap,S<:Space,T} <: AbstractMarkovMap
  m::M
  sp::S
  cache::Dict{Int,Array{T,2}}
end
function MarkovInverseCache(m::AbstractMarkovMap,s::Space)
  @assert rangedomain(m) == domain(s)
  MarkovInverseCache{typeof(m),typeof(s),eltype(domain(m))}(m,s,Dict{Int,Array{eltype(domain(m)),2}}())
end

MarkovInverseCache(m::AbstractMarkovMap) = MarkovInverseCache(m,ApproxFun.Space(rangedomain(m)))

ApproxFun.domain(mc::MarkovInverseCache) = domain(mc.m)
rangedomain(mc::MarkovInverseCache) = rangedomain(mc.m)
length(mc::MarkovInverseCache) = length(mc.m)
rangespace(mc::MarkovInverseCache) = mc.sp
function extenddata!(m::MarkovInverseCache,n::Int)
  v = Array(eltype(m),length(m),n)
  pts = points(rangespace(m),n)
  for k = 1:n
    for i = 1:length(m)
      v[i,k] = mapinv(m.m,i,pts[k])
    end
  end
  push!(m.cache,n=>v)
end

for FN in (:mapinv,:mapinvD)
  @eval $FN(m::MarkovInverseCache,i::Integer,x::Number) = $FN(m.m,i,x)
end

function mapinv(m::MarkovInverseCache,i::Integer,x::InterpolationNode)
  rangespace(m) != space(x) && return mapinv(m,i,convert(eltype(x),x))
  ~haskey(m.cache,x.n) && extenddata!(m,x.n)
  m.cache[x.n][i,x.k]
end
# mapinvD{M<:InverseDerivativeMarkovMap}(m::MarkovInverseCache{M},i::Integer,x::InterpolationNode) =
#   mapinvD(m.m,i,x)
mapinvD(m::MarkovInverseCache,i::Integer,x::InterpolationNode) =
  isa(m.branches[i],RevExpandingBranch) ? mapinvD(m.m,i,mapinv(m,i,x)) : 1/mapD(m.m,i,mapinv(m,i,x))
# mapDsign(m::MarkovInverseCache,i::Integer) = mapDsign(m.m,i)
