export NeutralBranch

# ExpandingBranch
@compat abstract type AbstractBranch{D<:Domain,R<:Domain}; end

@compat abstract type ExpandingBranch{D<:Domain,R<:Domain} <: AbstractBranch{D,R}; end

Base.summary(b::ExpandingBranch) =  string(typeof(b).name.name)*":"*string(domain(b))*"↦"*string(rangedomain(b)) #branches??
ApproxFun.prectype(b::ExpandingBranch) = prectype(rangedomain(b))
eltype(b::ExpandingBranch) = eltype(rangedomain(b))
# Base.show(io::IO,b::ExpandingBranch) = print(io,typeof(b)) #temporary

@compat struct FwdExpandingBranch{ff,gg,D<:Domain,R<:Domain} <: ExpandingBranch{D,R}
  f::ff
  dfdx::gg
  domain::D
  rangedomain::R
  #  sgn::T
  # function FwdExpandingBranch(fc,dfdxc,dom,ran)
    # @assert all([in(fc(p),∂(ran)) for p in ∂(dom)])
  #   new(fc,dfdxc,dom,ran)
  # end
end
function FwdExpandingBranch(f,dfdx,dom,ran)
  domd = convert(Domain,dom); randm = convert(Domain,ran)
  FwdExpandingBranch{typeof(f),typeof(dfdx),typeof(domd),typeof(randm)}(f,dfdx,domd,randm)
end

unsafe_call(b::FwdExpandingBranch,x) = b.f(x)
unsafe_mapD(b::FwdExpandingBranch,x) = b.dfdx(x)
unsafe_mapP(b::FwdExpandingBranch,x) = (unsafe_call(b,x),unsafe_mapD(b,x))
unsafe_mapinv(b::FwdExpandingBranch,x) = domain_newton(b.f,b.dfdx,x,b.domain,domain_guess(x,b.domain,b.rangedomain))
unsafe_mapinvD(b::FwdExpandingBranch,x) = inv(b.dfdx(unsafe_mapinv(b,x)))
function unsafe_mapinvP(b::FwdExpandingBranch,x)
  vx = unsafe_mapinv(b,x)
  (vx,inv(unsafe_mapD(b,vx)))
end


# RevExpandingBranch

@compat struct RevExpandingBranch{ff,gg,D<:Domain,R<:Domain} <: ExpandingBranch{D,R}
  v::ff
  dvdx::gg
  domain::D
  rangedomain::R
  #   sgn::T
  # function RevExpandingBranch(vc,dvdxc,dom,ran)
  #   new(vc,dvdxc,dom,ran)
  # end
end
function RevExpandingBranch(v,dvdx,dom,ran)
    domd = convert(Domain,dom); randm = convert(Domain,ran)
    RevExpandingBranch{typeof(v),typeof(dvdx),typeof(domd),typeof(randm)}(v,dvdx,domd,randm)
end


unsafe_call(b::RevExpandingBranch,x) = domain_newton(b.v,b.dvdx,x,b.rangedomain,domain_guess(x,b.rangedomain,b.domain))
unsafe_mapD(b::RevExpandingBranch,x) =  inv(b.dvdx(unsafe_call(b,x)))
function unsafe_mapP(b::RevExpandingBranch,x)
  fx = unsafe_call(b,x)
  (fx,inv(unsafe_mapinvD(b,fx)))
end
unsafe_mapinv(b::RevExpandingBranch,x) = b.v(x)
unsafe_mapinvD(b::RevExpandingBranch,x) = b.dvdx(x)
unsafe_mapinvP(b::RevExpandingBranch,x) = (unsafe_mapinv(b,x),unsafe_mapinvD(b,x))

for B in (FwdExpandingBranch,RevExpandingBranch)
  for (M,UNS_M) in ((:mapinv,:unsafe_mapinv),(:mapinvD,:unsafe_mapinvD),(:mapinvP,:unsafe_mapinvP))
    @eval $M(b::$B,x) = begin @assert in(x,rangedomain(b)); $UNS_M(b,x); end
  end

  # we are assuming that domains are unoriented
  @eval (b::$B)(x::Segment) = Segment(sort(b.(∂(x)))...)
  #@eval (b::$B)(x::Interval) = Interval(sort(b.(extrema(x)))...)
end

mapinv(b::ExpandingBranch,x::Segment) = Segment((mapinv.(b,∂(x)))...)
# mapinv(b::ExpandingBranch,x::Interval) = Interval(sort(mapinv.(b,extrema(x)))...)


# UNDER CONSTRUCTION: NeutralBranch
@compat struct NeutralBranch
end

# branch constructors

autodiff(f,d) = autodiff_dual(f,ApproxFun.checkpoints(convert(Domain,d)))
autodiff(f::Fun,d) = f'

function autodiff_dual(f,bi)
  fd = FunctionDerivative(f)
  try
    for b in bi
      fd(b)
    end
  catch e
    isa(e,MethodError) && error("To use automatic differentiation, your function must accept DualNumbers")
    throw(e)
  end
  fd
end

@compat const DomainInput = Union{Domain,IntervalSets.AbstractInterval}

function branch(f,dom,ran,diff=autodiff(f,(dir==Forward ? dom : ran));dir=Forward,
                      ftype=typeof(f),difftype=typeof(diff))
  domd = convert(Domain,dom); randm  = convert(Domain,ran);
  dir==Forward ? FwdExpandingBranch{ftype,difftype,typeof(domd),typeof(randm)}(f,diff,domd,randm) :
        RevExpandingBranch{ftype,difftype,typeof(domd),typeof(randm)}(f,diff,domd,randm)
end
branch(f,dom,ran,diff::Nothing; dir =Forward) = branch(f,dom,ran;dir=dir)

@deprecate branch(f,dfdx,dom::Domain,ran::Domain;dir=Forward) branch(f,dom,ran,diff;dir=dir)


for TYP in (:FwdExpandingBranch,:RevExpandingBranch,:NeutralBranch)
  @eval (b::$TYP)(x) = temp_in(x,b.domain) ? unsafe_call(b,x) : error("DomainError: $x ∉ $(b.domain)")
  @eval mapD(b::$TYP,x) = temp_in(x,b.domain) ? unsafe_mapD(b,x) : error("DomainError: $x ∉ $(b.domain)")
  @eval mapP(b::$TYP,x) = temp_in(x,b.domain) ? unsafe_mapP(b,x) : error("DomainError: $x ∉ $(b.domain)")
  @eval domain(b::$TYP) = b.domain
  @eval rangedomain(b::$TYP) = b.rangedomain
end

# transferbranch

function transferbranch_int_edges(x,y,b::ExpandingBranch)
  x ∈ rangedomain(b) && y ∈ rangedomain(b) && (return x,y)
  iv = Segment(x,y)∩rangedomain(b)
  iv.a, iv.b
end

function transferbranch(x,b::ExpandingBranch,f)
  x ∉ rangedomain(b) && return zero(typeof(x))
  (v,dvdx) = mapinvP(b,x)
  abs(det(dvdx))*f(v)
end
function transferbranch_int(x,y,b::ExpandingBranch,f)
  x ∉ rangedomain(b) && y ∉ rangedomain(b) && (return zero(typeof(x)))
  x,y = transferbranch_int_edges(x,y,b)
  csf = cumsum(f)
  vy = unsafe_mapinv(b,y); vx = unsafe_mapinv(b,x)
  sgn = sign((vy-vx)/(y-x))
  sgn*(csf(vy)-csf(vx))
end

function transferbranch(x,b::ExpandingBranch,sk::BasisFun)
  x ∉ rangedomain(b) && return zero(typeof(x))
  (v,dvdx) = unsafe_mapinvP(b,x)
  abs(det(dvdx))*getbasisfun(v,sk)
end
function transferbranch_int(x,y,b::ExpandingBranch,sk::BasisFun)
  x ∉ rangedomain(b) && y ∉ rangedomain(b) && (return zero(typeof(x)))
  x,y = transferbranch_int_edges(x,y,b)
  vy = unsafe_mapinv(b,y); unsafe_vx = unsafe_mapinv(b,x)
  sgn = sign((vy-vx)/(y-x))
  sgn*(getbasisfun_int(vy,sk)-getbasisfun_int(vx,sk))
end
