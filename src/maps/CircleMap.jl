# CircleMaps

@compat abstract type AbstractCircleMap{D<:PeriodicDomain,R<:PeriodicDomain} <: AbstractMarkovMap{D,R} end

@compat struct FwdCircleMap{D<:PeriodicDomain,R<:PeriodicDomain,ff,gg,T} <: AbstractCircleMap{D,R}
  f::ff
  dfdx::gg
  domain::D
  rangedomain::R
  cover::Int
  fa::T
  fb::T
end
function FwdCircleMap(f::ff,domd::D,randm::R,dfdx::gg=autodiff(f,domd)) where {D<:PeriodicDomain,R<:PeriodicDomain,ff,gg}
  # @assert isempty(∂(randm)) isempty(∂(domd))
  fa = f(first(domd)); fb = f(last(domd))
  cover_est = (fb-fa)/arclength(randm)
  cover_integer = round(Int,cover_est)
  cover_est ≈ cover_integer || error("Circle map lift does not have integer covering number.")
  FwdCircleMap{D,R,ff,gg,typeof(fa)}(f,dfdx,domd,randm,abs(cover_integer),fa,fb)
end
function FwdCircleMap(f::ff,dom,ran,dfdx::gg=autodiff(f,convert(PeriodicDomain,dom))) where {ff,gg}
  domd = convert(PeriodicDomain,dom); randm = convert(PeriodicDomain,ran)
  FwdCircleMap{typeof(domd),typeof(randm),ff,gg}(f,domd,randm,dfdx)
end

(m::FwdCircleMap)(x) = domain_mod(m.f(x),m.rangedomain)
mapD(m::FwdCircleMap,x) = m.dfdx(x)
mapP(m::FwdCircleMap,x) = (m(x),m.dfdx(x))

mapinv(m::FwdCircleMap,i::Integer,y) = domain_mod(domain_newton(m.f,m.dfdx,interval_mod(y+i*arclength(m.rangedomain),m.fa,m.fb),
    m.domain),m.domain)
mapinvD(m::FwdCircleMap,i::Integer,y) = m.dfdx(mapinv(m,i,y))
function mapinvP(m::FwdCircleMap,i::Integer,y)
  x = mapinv(m,i,y)
  (x,inv(m.dfdx(x)))
end


@compat struct RevCircleMap{D<:Domain,R<:Domain,ff,gg,T} <: AbstractCircleMap{D,R}
  v::ff
  dvdx::gg
  domain::D
  rangedomain::R
  cover::Int
  va::T
  vb::T
end

function RevCircleMap(v::ff,domd::D,randm::R,dvdx::gg=autodiff(v,randm), maxcover=10000) where {D<:Domain,R<:Domain,ff,gg}
  # @assert isempty(∂(ran)) isempty(∂(dom))
  ra = first(randm); va = v(ra)
  dr = arclength(randm); dd = arclength(domd)
  cover = 1; vr = v(ra+dr)-va
  while (abs(vr) <= dd && cover < maxcover)
    vr = v(ra+cover*dr)-va
    abs(vr) ≈ dd && break
    abs(vr) > dd && error("Inverse lift doesn't appear to have an inverse")
    cover += 1
  end
  cover == maxcover && error("Can't get to the end of the inverse lift after $maxcover steps")

  RevCircleMap{D,R,ff,gg,typeof(va)}(v,dvdx,domd,randm,cover,va,v(last(randm)))
end
function RevCircleMap(v::ff,dom,ran=dom,dvdx::gg=autodiff(v,ran)) where {ff,gg}
  domd = convert(PeriodicDomain,dom); randm = convert(PeriodicDomain,ran)
  RevCircleMap{typeof(domd),typeof(randm),ff,gg}(v,domd,randm,dvdx)
end

mapL(m::RevCircleMap,x) = domain_newton(m.v,m.dvdx,interval_mod(x,m.va,m.vb),m.rangedomain)
(m::RevCircleMap)(x) = domain_mod(mapL(m,x),m.domain)
mapD(m::RevCircleMap,x) = inv(m.dvdx(mapL(m,x)))
function mapP(m::RevCircleMap,x)
  y = mapL(m,x)
  (interval_mod(y,m.domain),inv(m.dvdx(y)))
end

mapinv(m::RevCircleMap,i::Integer,y) = m.v(y+i*arclength(m.rangedomain))
mapinvD(m::RevCircleMap,i::Integer,y) = m.dvdx(y+i*arclength(m.rangedomain))
function mapinvP(m::RevCircleMap,i::Integer,y)
  ys = y+i*arclength(m.rangedomain)
  (m.v(ys),m.dvdx(ys))
end




# autodiff
# FwdCircleMap(f,d,r=d;diff=autodiff(f,d)) = FwdCircleMap(f,diff,d,r);
# RevCircleMap(f,d,r=d;diff=autodiff(f,r)) = RevCircleMap(f,diff,d,r);
CircleMap(f,d,r=d;dir=Forward,diff=autodiff(f,dir==Forward ? d : r)) = dir == Forward ?
  FwdCircleMap(f,d,r,diff) : RevCircleMap(f,d,r,diff)

ncover(m::AbstractCircleMap) = m.cover
eachbranchindex(m::AbstractCircleMap) = 1:m.cover

for TYP in (:FwdCircleMap,:RevCircleMap)
  @eval ApproxFun.domain(m::$TYP) = m.domain
  @eval rangedomain(m::$TYP) = m.rangedomain
end

# Transfer function

function transferfunction(x,m::AbstractCircleMap,f)
  y = zero(eltype(x));
  for b = 1:ncover(m)
    (v,dvdx) = mapinvP(m,b,x)
    y += abs(det(dvdx))*f(v)
  end;
  y
end

function transferfunction_int(x,y,m::AbstractCircleMap,f)
  q = zero(eltype(x));
  csf = cumsum(f)
  for b = 1:ncover(m)
    vy = mapinv(m,b,y); vx = mapinv(m,b,x)
    sgn = sign((vy-vx)/(y-x))
    q += sgn*(csf(vy)-csf(vx))
  end;
  q
end

function transferfunction(x,m::AbstractCircleMap,sk::BasisFun)
  y = zero(eltype(x));
  for b = 1:ncover(m)
    (v,dvdx) = mapinvP(m,b,x)
    y += abs(det(dvdx)).*getbasisfun(v,sk)
  end
  y
end
function transferfunction_int(x,y,m::AbstractCircleMap,sk::BasisFun)
  q = zero(eltype(x));
  for b = 1:ncover(m)
    vy = mapinv(m,b,y); vx = mapinv(m,b,x)
    sgn = sign((vy-vx)/(y-x))
    q += sgn*(getbasisfun_int(vy,sk)-getbasisfun_int(vx,sk))
  end
  q
end
