# CircleMaps

abstract AbstractCircleMap{D<:Domain,R<:Domain} <: AbstractMarkovMap{D,R}

immutable FwdCircleMap{D<:Domain,R<:Domain,ff,gg,T} <: AbstractCircleMap{D,R}
  f::ff
  dfdx::gg
  domain::D
  rangedomain::R
  cover::Int
  fa::T
  fb::T

  function FwdCircleMap(f,domd,randm,dfdx)
    # @assert isempty(∂(randm)) isempty(∂(domd))
    fa = f(first(domd)); fb = f(last(domd))
    cover_est = (fb-fa)/arclength(randm)
    cover_integer = round(Int,cover_est)
    cover_est ≈ cover_integer || error("Circle map lift does not have integer covering number.")
    new(f,dfdx,domd,randm,abs(cover_integer),fa,fb)
  end
end
function FwdCircleMap{ff,gg}(f::ff,dom,ran,dfdx::gg=autodiff(f,dom))
  domd = PeriodicDomain(dom); randm = PeriodicDomain(ran)
  FwdCircleMap{typeof(domd),typeof(randm),ff,gg,eltype(domd)}(f,domd,randm,dfdx)
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


immutable RevCircleMap{D<:Domain,R<:Domain,ff,gg,T} <: AbstractCircleMap{D,R}
  v::ff
  dvdx::gg
  domain::D
  rangedomain::R
  cover::Int
  va::T
  vb::T

  function RevCircleMap(v,domd,randm,dvdx)

    # @assert isempty(∂(ran)) isempty(∂(dom))
    ra = first(randm); va = v(ra)
    dr = arclength(randm); dd = arclength(domd)
    cover = 1; vr = v(ra+dr)-va
    while (abs(vr) <= dd && cover < 10000)
      vr = v(ra+cover*dr)-va
      abs(vr) ≈ dd && break
      abs(vr) > dd && error("Inverse lift doesn't appear to have an inverse")
      cover += 1
    end
    cover == 10000 && error("Can't get to the end of the inverse lift after 10000 steps")

    new(v,dvdx,domd,randm,cover,va,v(last(randm)))
  end
end
function RevCircleMap{ff,gg}(v::ff,dom,ran=dom,dvdx::gg=autodiff(v,ran))
  domd = PeriodicDomain(dom); randm = PeriodicDomain(ran)
  RevCircleMap{typeof(domd),typeof(randm),ff,gg,eltype(randm)}(v,domd,randm,dvdx)
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
CircleMap(f,d,r=d;dir=Forward,diff=autodiff(f,dir=Forward ? d : r)) = dir == Forward ?
  FwdCircleMap(f,d,r,diff) : RevCircleMap(f,d,r,diff)

ncover(m::AbstractCircleMap) = m.cover

for TYP in (:FwdCircleMap,:RevCircleMap)
  @eval ApproxFun.domain(m::$TYP) = m.domain
  @eval rangedomain(m::$TYP) = m.rangedomain
end
