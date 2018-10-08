struct ComposedMarkovMap{T<:Tuple,D<:Domain,R<:Domain} <: AbstractMarkovMap{D,R}
  maps::T
  domain::D
  rangedomain::R
end
(∘)(f::AbstractMarkovMap,g::AbstractMarkovMap) = ComposedMarkovMap(f,g)

function ComposedMarkovMap(maps...)
  ran = rangedomain(maps[1])
  mps = isa(maps[1],ComposedMarkovMap) ? maps[1].maps : (maps[1],)
  for i = 1:length(maps)-1
    @assert rangedomain(maps[i+1]) == domain(maps[i])
    mps = isa(maps[i+1],ComposedMarkovMap) ? (mps...,maps[i+1].maps...) : (mps...,maps[i+1])
  end
  dom = domain(maps[end])

  ComposedMarkovMap(mps,dom,ran)
end

complength(c::ComposedMarkovMap) = length(c.maps)

function (c::ComposedMarkovMap)(x)
  f = c.maps[end](x)
  for i = complength(c)-1:-1:1
    f = c.maps[i](f)
  end
  f
end

function mapP(c::ComposedMarkovMap,x)
  f,dfdx = mapP(c.maps[end],x)
  for i = complength(c)-1:-1:1
    f,dfdxs = mapP(c.maps[i],f)
    dfdx *= dfdxs
  end
  f,dfdx
end
mapD(c::ComposedMarkovMap,x) = mapP(c,x)[2]

function mapinv(c::ComposedMarkovMap,b,x)
  v = mapinv(c.maps[1],b[1],x)
  for i = 2:complength(c)
    v = mapinv(c.maps[i],b[i],v)
  end
  v
end
function mapinvP(c::ComposedMarkovMap,b,x)
  v,dvdx = mapinvP(c.maps[1],b[1],x)
  for i = 2:complength(c)
    v,dvdxs = mapinvP(c.maps[i],b[i],v)
    dvdx *= dvdxs
  end
  v,dvdx
end
mapinvD(c::ComposedMarkovMap,b,x) = mapinvP(c,b,x)[2]

# TODO: map(P,D)(c,b,x)

function getbranchind(m::ComposedMarkovMap,x)
  temp_in(x,m.domain) || error("DomainError: $x ∉ $(m.domain)")
  fx = x
  br = getbranchind(m.maps[end],fx)
  for i = complength(c)-1:-1:1
    fx = m.maps[i+1](fx)
    br = (getbranchind(m.maps[i],fx),br...)
  end
  br
end

nbranches(C::ComposedMarkovMap) = prod(nbranches(mm for mm in C.maps))
eachbranchindex(C::ComposedMarkovMap) = product(eachbranchindex(mm) for mm in C.maps)
branchindtype(C::ComposedMarkovMap) = Array{promote_type([branchindtype(mm) for mm in C.maps]...)}
nneutral(C::ComposedMarkovMap) = 0 # assume we can't do this
neutralfixedpoints(C::ComposedMarkovMap) = []

#TODO: must be faster??
<<<<<<< HEAD
function transferfunction(x,m::ComposedMarkovMap{Tuple{M}},f,T) where {M<:AbstractMarkovMap}
  transferfunction(x,m.maps[1],f,T)
=======
function transferfunction{M<:AbstractMarkovMap}(x,m::ComposedMarkovMap{Tuple{M}},f)
  transferfunction(x,m.maps[1],f)
>>>>>>> 4c33d1554aa02b6a20399f0090f695de8ad3e517
end

struct TransferCall{M<:AbstractMarkovMap,ff}
  m::M
  f::ff
end
# TransferCall(m,f) = TransferCall{typeof(m),typeof(f),T}(m,f)
(t::TransferCall)(x) = transferfunction(x,t.m,t.f)

function transferfunction(x,m::ComposedMarkovMap,f)
  m2 = ComposedMarkovMap(m.maps[2:end],domain(m.maps[end]),rangedomain(m.maps[2]))
  transferfunction(x,m.maps[1],#x->transferfunction(x,m2,f))
    TransferCall(m2,f))
end

# TODO: transferfunction_int


# TODO
  struct ComposedCircleMap{M<:AbstractMarkovMap,D<:Domain,R<:Domain} <: AbstractMarkovMap{D,R}
    maps::NTuple{M}
    domain::D
    rangedomain::R
    end


    # stuff you can do for both
ComposedMap{M,D<:Domain,R<:Domain} = Union{ComposedMarkovMap{M,D,R},ComposedCircleMap{M,D,R}}

for FUN in (:domain,:rangedomain)
  @eval ($FUN)(c::ComposedMarkovMap) = c.$FUN
end
