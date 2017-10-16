export NeutralBranch

# MarkovBranch

@compat abstract type MarkovBranch{D<:Domain,R<:Domain}; end

Base.summary(b::MarkovBranch) =  string(typeof(b).name.name)*":"*string(domain(b))*"↦"*string(rangedomain(b)) #branches??
Base.eltype(b::MarkovBranch) = eltype(rangedomain(b))
# Base.show(io::IO,b::MarkovBranch) = print(io,typeof(b)) #temporary

@compat struct FwdExpandingBranch{ff,gg,D<:Domain,R<:Domain} <: MarkovBranch{D,R}
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
function FwdExpandingBranch{D,R,ff,gg}(f::ff,dfdx::gg,dom::D,ran::R)
  domd = Domain(dom); randm = Domain(ran)
  FwdExpandingBranch{typeof(f),typeof(dfdx),typeof(domd),typeof(randm)}(f,dfdx,domd,randm)
end

unsafe_call(b::FwdExpandingBranch,x) = b.f(x)
unsafe_mapD(b::FwdExpandingBranch,x) = b.dfdx(x)
unsafe_mapP(b::FwdExpandingBranch,x) = (unsafe_call(b,x),unsafe_mapD(b,x))
mapinv(b::FwdExpandingBranch,x) = domain_newton(b.f,b.dfdx,x,b.domain,domain_guess(x,b.domain,b.rangedomain))
mapinvD(b::FwdExpandingBranch,x) = inv(b.dfdx(mapinv(b,x)))
function mapinvP(b::FwdExpandingBranch,x)
  vx = mapinv(b,x)
  (vx,inv(unsafe_mapD(b,vx)))
end


# RevExpandingBranch

@compat struct RevExpandingBranch{ff,gg,D<:Domain,R<:Domain} <: MarkovBranch{D,R}
  v::ff
  dvdx::gg
  domain::D
  rangedomain::R
  #   sgn::T
  # function RevExpandingBranch(vc,dvdxc,dom,ran)
  #   new(vc,dvdxc,dom,ran)
  # end
end
function RevExpandingBranch{D,R,ff,gg}(v::ff,dvdx::gg,dom::D,ran::R)
    domd = Domain(dom); randm = Domain(ran)
    RevExpandingBranch{typeof(v),typeof(dvdx),typeof(domd),typeof(randm)}(v,dvdx,domd,randm)
end


unsafe_call(b::RevExpandingBranch,x) = domain_newton(b.v,b.dvdx,x,b.rangedomain,domain_guess(x,b.rangedomain,b.domain))
unsafe_mapD(b::RevExpandingBranch,x) =  inv(b.dvdx(unsafe_call(b,x)))
function unsafe_mapP(b::RevExpandingBranch,x)
  fx = unsafe_call(b,x)
  (fx,inv(mapinvD(b,fx)))
end
mapinv(b::RevExpandingBranch,x) = b.v(x)
mapinvD(b::RevExpandingBranch,x) = b.dvdx(x)
mapinvP(b::RevExpandingBranch,x) = (mapinv(b,x),mapinvD(b,x))



# UNDE CONSTRUCTION: NeutralBranch
@compat struct NeutralBranch
end

# branch constructors

autodiff(f,d) = autodiff_dual(f,ApproxFun.checkpoints(Domain(d)))
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

function branch(f,dom,ran,diff=autodiff(f,(dir=Forward ? dom : ran));dir=Forward,
                      ftype=typeof(f),difftype=typeof(diff))
  domd = Domain(dom); randm  = Domain(ran);
  dir==Forward ? FwdExpandingBranch{ftype,difftype,typeof(domd),typeof(randm)}(f,diff,domd,randm) :
        RevExpandingBranch{ftype,difftype,typeof(domd),typeof(randm)}(f,diff,domd,randm)
end
branch(f,dom,ran,diff::Void;dir=Forward) = branch(f,dom,ran;dir)

@deprecate branch(f,dfdx,dom::Domain,ran::Domain;dir=Forward) branch(f,dom,ran,diff;dir=dir)


for TYP in (:FwdExpandingBranch,:RevExpandingBranch,:NeutralBranch)
  @eval (b::$TYP)(x) = temp_in(x,b.domain) ? unsafe_call(b,x) : error("DomainError: $x ∉ $(b.domain)")
  @eval mapD(b::$TYP,x) = temp_in(x,b.domain) ? unsafe_mapD(b,x) : error("DomainError: $x ∉ $(b.domain)")
  @eval mapP(b::$TYP,x) = temp_in(x,b.domain) ? unsafe_mapP(b,x) : error("DomainError: $x ∉ $(b.domain)")
  @eval domain(b::$TYP) = b.domain
  @eval rangedomain(b::$TYP) = b.rangedomain
end

# transferbranch

function transferbranch(x,b::MarkovBranch,f,T)
  (v,dvdx) = mapinvP(b,x)
  abs(det(dvdx))*f(v)
end
function transferbranch_int(x,y,b::MarkovBranch,f,T)
  csf = cumsum(f)
  vy = mapinv(b,y); vx = mapinv(b,x)
  sgn = sign((vy-vx)/(y-x))
  sgn*(csf(vy)-csf(vx))
end

function transferbranch(x,b::MarkovBranch,sk::BasisFun,T)
  (v,dvdx) = mapinvP(b,x)
  abs(det(dvdx))*getbasisfun(v,sk,T)
end
function transferbranch_int(x,y,b::MarkovBranch,sk::BasisFun,T)
  vy = mapinv(b,y); vx = mapinv(b,x)
  sgn = sign((vy-vx)/(y-x))
  sgn*(getbasisfun_int(vy,sk,T)-getbasisfun_int(vx,sk,T))
end
