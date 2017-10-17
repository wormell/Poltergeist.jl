struct ComposedMarkovMap{M<:AbstractMarkovMap,D<:Domain,R<:Domain,N} <: AbstractMarkovMap{D,R}
  maps::NTuple{N,M}
  domain::D
  rangedomain::R
end

function ComposedMarkovMap(maps...)
  T = typejoin(map(typeof,maps)...)
  N = length(maps)
  ran = rangedomain(maps[1])
  for i = 1:N-1
    @assert rangedomain(maps[i+1]) == domain(maps[i])
  end
  dom = domain(maps[end])
  ComposedMarkovMap{T,typeof(dom),typeof(ran),N}(convert(NTuple{N,T},maps),dom,ran)
end
(∘)(f::AbstractMarkovMap,g::AbstractMarkovMap) = ComposedMarkovMap(f,g)

complength{M,D,R,N}(::ComposedMarkovMap{M,D,R,N}) = N

function (c::ComposedMarkovMap{M,D,R,N})(x) where {M,D,R,N}
  f = c.maps[end](x)
  for i = N-1:-1:1
    f = c.maps[i](f)
  end
  f
end

function mapP(c::ComposedMarkovMap{M,D,R,N},x) where {M,D,R,N}
  f,dfdx = mapP(c.maps[end],x)
  for i = N-1:-1:1
    f,dfdxs = mapP(c.maps[i],f)
    dfdx *= dfdxs
  end
  f,dfdx
end
mapD(c::ComposedMarkovMap,x) = mapP(c,x)[2]

function mapinv(c::ComposedMarkovMap{M,D,R,N},b,x) where {M,D,R,N}
  v = mapinv(c.maps[1],b[1],x)
  for i = 2:N
    v = mapinv(c.maps[i],b[i],v)
  end
  v
end
function mapinvP(c::ComposedMarkovMap{M,D,R,N},b,x) where {M,D,R,N}
  v,dvdx = mapinvP(c.maps[1],b[1],x)
  for i = 2:N
    v,dvdxs = mapinvP(c.maps[i],b[i],v)
    dvdx *= dvdxs
  end
  v,dvdx
end
mapinvD(c::ComposedMarkovMap,b,x) = mapinvP(c,b,x)[2]

# TODO: map(P,D)(c,b,x)

function getbranch(m::ComposedMarkovMap{M,D,R,N},x) where {M,D,R,N}
  temp_in(x,m.domain) || error("DomainError: $x ∉ $(m.domain)")
  fx = x
  br = getbranch(m.maps[N],fx)
  for i = N-1:-1:1
    fx = m.maps[i+1](fx)
    br = (getbranch(m.maps[i],fx),br...)
  end
  br
end

nbranches(C::ComposedMarkovMap) = prod(nbranches(mm for mm in C.maps))

function transferfunction{M,D,R,N}(x,m::ComposedMarkovMap{M,D,R,N},f,T)
  m2 = ComposedMarkovMap(m.maps[2:end],m.domain,rangedomain(m.maps[2]))
  transferfunction(x,m.maps[1],x->transferfunction(x,m2,f,T),T)
end
function transferfunction{M,D,R}(x,m::ComposedMarkovMap{M,D,R,1},f,T)
  transferfunction(x,m.maps[1],f,T)
end


# TODO: transferbranch_int


# TODO
  struct ComposedCircleMap{M<:AbstractMarkovMap,D<:Domain,R<:Domain} <: AbstractMarkovMap{D,R}
    maps::NTuple{M}
    domain::D
    rangedomain::R
    end


    # stuff you can do for both
ComposedMap{M,D<:Domain,R<:Range} = Union{ComposedMarkovMap{M,D,R},ComposedCircleMap{M,D,R}}

for FUN in (:domain,:rangedomain)
  @eval ($FUN)(c::ComposedMarkovMap) = c.$FUN
end
