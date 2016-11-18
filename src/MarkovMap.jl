export domain, rangedomain

abstract AbstractMarkovMap{D<:Domain,R<:Domain,T,FF} <: Function
#abstract AbstractDerivativeMarkovMap{D<:Domain,R<:Domain,T,FF} <: AbstractMarkovMap{D,R,T,FF}

Base.summary(m::AbstractMarkovMap) =  string(typeof(m).name.name)*":"*string(domain(m))*"↦"*string(rangedomain(m)) #branches??
Base.eltype(m::AbstractMarkovMap) = eltype(rangedomain(m))
Base.show(io::IO,m::AbstractMarkovMap) = print(io,typeof(m)) #temporary


immutable MarkovMap{D<:Domain,R<:Domain,T,ff} <: AbstractMarkovMap{D,R,T,ff}
  domain::D
  rangedomain::R
  f::Vector{ff}
  dfdx::Vector{ff}
  bl::Vector{T} #lower edges of branches
  bu::Vector{T} #upper edges of branches
  sgns::Vector{T} #sign of derivative of f on branches
  #     λhat::T #maximum derivative
  #     C1::T #first distortion bound
  # cache of inverse branches of f

end
# MarkovMap{T,ff<:Number}(dom::Domain,ran::Domain,f::Vector{ff},dfdx::Vector{ff},bl::Vector{T},bu::Vector{T}) =
#   MarkovMap{T,ff}(dom,ran,f,dfdx,bl,bu)

function MarkovMap{D,R,T,ff}(
    dom::D,ran::R,f::Vector{ff},dfdx::Vector{ff},bl::Vector{T},bu::Vector{T})
  @assert length(bl) == length(bu)
  @assert length(f) == length(bl)
  @assert length(dfdx) == length(bl)
  #        @assert length(bl) == length(sgns)
  #        @assert all(sgns.^2 == 1)
  @assert all([in(bll,dom) for bll in bl])
  @assert all([in(buu,dom) for buu in bu])
  @assert all(bl.<bu)
  bl_sortperm = sortperm(bl)
  @assert all(bu[bl_sortperm[1:end-1]] .<= bl[bl_sortperm[2:end]])
  #     if isa(ran,ApproxFun.IntervalDomain)
  #       @assert all([in(f[i](bl[i]),∂(ran)) for i = 1:length(bl)])
  #       @assert all([in(f[i](bu[i]),∂(ran)) for i = 1:length(bu)])
  #     end
  sgns = T[sign(f[i](bu[i])-f[i](bl[i])) for i = 1:length(bl)]
  MarkovMap(dom,ran,f,dfdx,bl,bu,sgns)
end

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

length(m::MarkovMap) = length(m.bl)
branch(m::MarkovMap,x::Number) = findfirst([x .<= bu[i] && x .>= bl[i] for i=1:length(bu)])
(m::MarkovMap)(i::Integer,x) = m.f[i](x)
(m::MarkovMap)(x::Number) = m.f[branch(m,x)](x)
markovD(m::MarkovMap,i::Integer,x) = m.dfdx[i](x)
markovD(m::MarkovMap,x::Number) = m.dfdx[branch(m,x)](x)
markovinv(m::MarkovMap,i::Integer,x::Number) = interval_newton(m.f[i],m.dfdx[i],x,m.bl[i],m.bu[i])
markovinvD(m::MarkovMap,i::Integer,x::Number) = 1/m.dfdx[i](markovinv(m,i,x))
markovDsign(m::MarkovMap,i::Integer) = m.sgns[i]

immutable InverseDerivativeMarkovMap{D<:Domain,R<:Domain,T,FF} <: AbstractMarkovMap{D,R,T,FF}
  domain::D
  rangedomain::R
  v::Vector{FF} #branches of inverse
  dvdx::Vector{FF} #derivatives of branches of inverse
  sgns::Vector{T} #sign of derivative of f on branches

  # assert stuff??
  #   function InverseDerivativeMarkovMap{FF,T}(dom::Domain,ran::Domain,v::Vector{FF},dvdx::Vector{FF},sgns::Vector{T})
  #       @assert length(v) == length(dvdx)
  #       @assert length(v) == length(sgns)
  #       @assert all(sgns.^2 == 1)
  #       new(dom,ran,v,dvdx,sgns)
  #   end
end

function InverseDerivativeMarkovMap{FF}(dom::Domain,ran::Domain,v::Vector{FF},dvdx::Vector{FF})
  sgnpt = rand(ran); #hack
  return InverseDerivativeMarkovMap(dom,ran,v,dvdx,eltype(dom)[sign(dvdxi(sgnpt)) for dvdxi in dvdx])
end
InverseDerivativeMarkovMap{FF}(dom::Domain,v::Vector{FF},dvdx::Vector{FF}) =
  InverseDerivativeMarkovMap(dom,dom,v,dvdx)


length(m::InverseDerivativeMarkovMap) = length(m.v)
markovinv(m::InverseDerivativeMarkovMap,i::Integer,x) = m.v[i](x)
markovinvD(m::InverseDerivativeMarkovMap,i::Integer,x) = m.dvdx[i](x)
markovDsign(m::InverseDerivativeMarkovMap,i::Integer) = m.sgns[i]

#getindex(m::InverseDerivativeMarkovMap,:dvdx) = m.dvdx

for TYP in (:MarkovMap,:InverseDerivativeMarkovMap)
  @eval ApproxFun.domain(m::$TYP) = m.domain
  @eval rangedomain(m::$TYP) = m.rangedomain
end

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
      v[i,k] = markovinv(m.m,i,pts[k])
    end
  end
  push!(m.cache,n=>v)
end

for FN in (:markovinv,:markovinvD)
  @eval $FN(m::MarkovInverseCache,i::Integer,x::Number) = $FN(m.m,i,x)
end

function markovinv(m::MarkovInverseCache,i::Integer,x::InterpolationNode)
  rangespace(m) != space(x) && return markovinv(m,i,convert(eltype(x),x))
  ~haskey(m.cache,x.n) && extenddata!(m,x.n)
  m.cache[x.n][i,x.k]
end
markovinvD{M<:InverseDerivativeMarkovMap}(m::MarkovInverseCache{M},i::Integer,x::InterpolationNode) =
  markovinvD(m.m,i,x)
markovinvD{M<:MarkovMap}(m::MarkovInverseCache{M},i::Integer,x::InterpolationNode) =
  1/markovD(m.m,i,markovinv(m,i,x))
markovDsign(m::MarkovInverseCache,i::Integer) = markovDsign(m.m,i)


# abstract MarkovRegularityWrapper{D<:Domain,R<:Domain,T,FF} <: AbstractMarkovMap{D,R,T,FF}

# type AnalyticMarkovRegularityWrapper{D<:Domain,R<:Domain,T,FF} <: MarkovRegularityWrapper{D,R,T,FF}
#   mm::AbstractMarkovMap{D,R,T,FF}
#   delta::T
#   p::T
#   K::T
# end

# function getanalyticp{D<:Domain,R<:PeriodicInterval,T,FF}(mm::AbstractDerivativeMarkovMap{D,R,T,FF},delta::T)
#   try
#     maxabs(Fun(x->imag(mm[:dvdx](x+im*delta)),space(R)))
#   catch
#     mmf = Fun(mm[:dvdx],space(R))
#     maxabs(Fun(x->imag(mmf(x+im*delta)),space(R)))
#   end
# end


# type DifferentiableMarkovRegularityWrapper{D<:Domain,R<:Domain,T,FF} <: MarkovRegularityWrapper{D,R,T,FF}
#   mm::AbstractMarkovMap{D,R,T,FF}
#   r::Int
#   W::Array{T}
# end

