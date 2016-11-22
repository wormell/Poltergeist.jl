export MarkovMap, branch, nbranches, induce

# MarkovBranch

abstract MarkovBranch{D<:Domain,R<:Domain}

Base.summary(b::MarkovBranch) =  string(typeof(b).name.name)*":"*string(domain(b))*"↦"*string(rangedomain(b)) #branches??
Base.eltype(b::MarkovBranch) = eltype(rangedomain(b))
Base.show(io::IO,b::MarkovBranch) = print(io,typeof(b)) #temporary

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
unsafe_mapP(b::FwdExpandingBranch,x::Number) = (unsafe_call(b,x),unsafe_mapD(b,x))
mapinv(b::FwdExpandingBranch,x::Number) = interval_newton(b.f,b.dfdx,x,b.domain,interval_guess(x,b.domain,b.rangedomain))
mapinvD(b::FwdExpandingBranch,x::Number) = 1/b.dfdx(mapinv(b,x))
function mapinvP(b::FwdExpandingBranch,x::Number)
  vx = mapinv(b,x)
  (vx,1/unsafe_mapD(b,vx))
end


# RevExpandingBranch

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

unsafe_call(b::RevExpandingBranch,x::Number) = interval_newton(b.v,b.dvdx,x,b.rangedomain,interval_guess(x,b.rangedomain,b.domain))
unsafe_mapD(b::RevExpandingBranch,x::Number) =  1/b.dvdx(unsafe_call(b,x))
function unsafe_mapP(b::RevExpandingBranch,x::Number)
  fx = unsafe_call(b,x)
  (fx,1/mapinvD(b,fx))
end
mapinv(b::RevExpandingBranch,x::Number) = b.v(x)
mapinvD(b::RevExpandingBranch,x::Number) = b.dvdx(x)
mapinvP(b::RevExpandingBranch,x::Number) = (mapinv(b,x),mapinvD(b,x))



# NeutralBranch

immutable NeutralBranch{A<:AbelFunction,D,R} <: MarkovBranch{D,R}
  ab::A
  domain::D
  rangedomain::R
  function NeutralBranch(abel,domain,rangedomain)
    @assert in(abel.p,∂(domain)∩∂(rangedomain))
    @assert issubset(domain,rangedomain)
    @assert domain != rangedomain
    @assert length(setdiff(∂(domain),∂(rangedomain)))==1
    @assert (setdiff(∂(domain),∂(rangedomain))[1]-abel.p)*abel.sgn > 0 "Abel function pointing in the wrong direction"
    @assert arclength(domain)*abel.h(arclength(domain)^abel.alpha) == arclength(rangedomain) "Branch not Markov"
    new(abel,domain,rangedomain)
  end
end
NeutralBranch{A<:AbelFunction,D<:Domain,R<:Domain}(abel::A,dom::D,ran::R) =
  NeutralBranch{typeof(abel),typeof(dom),typeof(ran)}(abel,dom,ran)


tocanonical(R::AbelFunction,x) = R.sgn*(x-R.p)
fromcanonical(R::AbelFunction,x) = R.p + R.sgn*x
for FUN in (:fromcanonical,:tocanonical)
  @eval $FUN(b::NeutralBranch,x) = $FUN(b.ab,x)
end

function NeutralBranch{hh,jj}(h::hh,dh::jj,alpha::Number,r0::Number,dom::Domain,ran::Domain)
  length(∂(dom)∩∂(ran)) != 1 && error("Edges of domain and range not compatible")
  p = (∂(dom)∩∂(ran))[1]
  sgn = sign(setdiff(∂(dom),∂(ran))[1]-p)
  abel = AbelFunction(h,dh,alpha,r0,p,sgn)
  zeroat!(abel,setdiff(∂(ran),∂(dom))[1])
  NeutralBranch(abel,dom,ran)
end

function unsafe_call(b::NeutralBranch,x)
  xc = tocanonical(b.ab,x)
  fromcanonical(b.ab,xc.*b.ab.h(xc.^b.ab.alpha))
end
function unsafe_mapD(b::NeutralBranch,x)
  xc = tocanonical(b.ab,x)
  alpha = b.ab.alpha
  b.ab.sgn*(b.ab.h(xc.^alpha) + alpha * xc.^alpha .* b.ab.dh(xc.^alpha))
end
unsafe_mapP(b::NeutralBranch,x) = (b(x),mapD(b,x))

function mapinv{T}(b::NeutralBranch,x::T)
  y = ringhconv(b.ab,toring(b.ab,x))::T
  fromring(b.ab,ringhconv(b.ab,hdisc_newton(b.ab,y)))
end
mapinvD(b::NeutralBranch,x) = mapD(b,mapinvD(b,x))
function mapinvP(b::NeutralBranch,x)
  vx = mapinv(b,x)
  (vx,1/unsafe_mapD(b,vx))
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

function autodiff_dual(f,bi::Number)
  fd = FunctionDerivative(f)
  try
    fd(bi)
  catch e
    isa(e,MethodError) && error("To use automatic differentiation, your function must accept DualNumbers")
    throw(e)
  end
  fd
end
autodiff_dual(f,bi::Domain) = autodiff_dual(f,rand(bi))

function branch{ff,T<:Union{Number,Domain}}(fs::Vector{ff},b::Vector{T},args...)
#   for (i,f) in enumerate(fs)
#     check_dualmethods(f,b[i])
#   end
  branch(fs,[autodiff_dual(fs[i],b[i]) for i in eachindex(fs)],b,args...)
end
branch{T<:Union{Number,Domain}}(f,b::T,args...) = branch(f,autodiff_dual(f,b),b,args...)

branch{ff<:Fun,T<:Union{Number,Domain}}(fs::Vector{ff},b::Vector{T},args...) =
  branch(fs,[f' for f in fs],b,args...)
branch{T<:Union{Number,Domain}}(f::Fun,b::T,args...) = branch(f,f',b,args...)




for TYP in (:FwdExpandingBranch,:RevExpandingBranch,:NeutralBranch)
  @eval (b::$TYP)(x::Number) = in(x,b.domain) ? unsafe_call(b,x) : throw(DomainError)
  @eval mapD(b::$TYP,x::Number) = in(x,b.domain) ? unsafe_mapD(b,x) : error("DomainError: $x in $(b.domain)")
  @eval mapP(b::$TYP,x::Number) = in(x,b.domain) ? unsafe_mapP(b,x) : throw(DomainError)
  @eval domain(b::$TYP) = b.domain
  @eval rangedomain(b::$TYP) = b.rangedomain
end




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

type InducedMarkovMap{M<:MarkovMap,B<:MarkovBranch,D<:Domain,R<:Domain} <: AbstractMarkovMap{D,R}
  domain::D
  rangedomain::R
  m::M
  b::B
  function InducedMarkovMap(dom,ran,m,b) #requires m,b: their domains -> m ∪ b
    @assert dom == ran
    @assert domain(m) == dom
    @assert rangedomain(b) == rangedomain(m)
    @assert domain(m) == setdiff(rangedomain(m),domain(b))
    new(dom,ran,m,b)
  end
end
InducedMarkovMap(dom::Domain,ran::Domain,m::MarkovMap,b::MarkovBranch) =
  InducedMarkovMap{typeof(m),typeof(b),typeof(dom),typeof(ran)}(dom,ran,m,b)


for TYP in (:MarkovMap,:InducedMarkovMap)
  @eval ApproxFun.domain(m::$TYP) = m.domain
  @eval rangedomain(m::$TYP) = m.rangedomain
end


# TODO: maybe roll into branch constructors? maybe remove??
function MarkovMap{D,R,ff}(
    dom::D,ran::R,f::Vector{ff},dfdx::Vector{ff})
  T = eltype(R)
  bl = Array(T,length(f)); bu = Array(T,length(f))
  tol = 10eps(maxabs([dom.a,dom.b]))
  for i = 1:length(f)
    ba = interval_newton(f[i],dfdx[i],ran.a,dom.a,dom.b,interval_guess(x,dom,ran),tol)
    bb = interval_newton(f[i],dfdx[i],ran.b,dom.a,dom.b,interval_guess(x,dom,ran),tol)
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
getbranch(m::MarkovMap,x::Number) = in(x,domain(m)) ? findfirst([in(x,domain(b)) for b in m.branches]) : throw(DomainError)

(m::MarkovMap)(i::Integer,x) = (m.branches[i])(x)
(m::MarkovMap)(x::Number) = (m.branches[getbranch(m,x)])(x)
for FUN in (:mapD,:mapP)
 @eval $FUN(m::MarkovMap,x::Number) = $FUN(m.branches[getbranch(m,x)],x)
end
for FUN in (:mapD,:mapP,:mapinv,:mapinvD,:mapinvP)
  @eval $FUN(m::MarkovMap,i::Integer,x::Number) = $FUN(m.branches[i],x)
end
#mapDsign(m::MarkovMap,i::Integer) = m.sgns[i]

#getindex(m::InverseDerivativeMarkovMap,:dvdx) = m.dvdx






# Inducing

function mapinduceP(b::MarkovBranch,x::Number) #hack
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
  InducedMarkovMap(dom,dom,MarkovMap(dom,m.rangedomain,m.branches[setdiff(1:nbranches(m),n)]),b)
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




# InverseCache

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

