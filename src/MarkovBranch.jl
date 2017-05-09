export NeutralBranch

# MarkovBranch

abstract MarkovBranch{D<:Domain,R<:Domain}

Base.summary(b::MarkovBranch) =  string(typeof(b).name.name)*":"*string(domain(b))*"↦"*string(rangedomain(b)) #branches??
Base.eltype(b::MarkovBranch) = eltype(rangedomain(b))
# Base.show(io::IO,b::MarkovBranch) = print(io,typeof(b)) #temporary

immutable FwdExpandingBranch{ff,gg,D<:Domain,R<:Domain} <: MarkovBranch{D,R}
  f::ff
  dfdx::gg
  domain::D
  rangedomain::R
  #  sgn::T
  function FwdExpandingBranch(fc,dfdxc,dom,ran)
    # @assert all([in(fc(p),∂(ran)) for p in ∂(dom)])
    new(fc,dfdxc,dom,ran)
  end
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
    @assert arclength(domain)*abel.h(arclength(domain)^abel.alpha) ≈ arclength(rangedomain) "Branch not Markov"
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
  (vx,inv(unsafe_mapD(b,vx)))
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

typealias DomainInput Union{Domain,IntervalSets.AbstractInterval}

function branch(f,dom,ran,diff=autodiff(f,(dir=Forward ? dom : ran));dir=Forward,
                      ftype=typeof(f),difftype=typeof(diff))
  domd = Domain(dom); randm  = Domain(ran);
  dir==Forward ? FwdExpandingBranch{ftype,difftype,typeof(domd),typeof(randm)}(f,diff,domd,randm) :
        RevExpandingBranch{ftype,difftype,typeof(domd),typeof(randm)}(f,diff,domd,randm)
end
branch(f,dom,ran,diff::Void;dir=Forward) = branch(f,dom,ran;dir)

@deprecate branch(f,dfdx,dom::Domain,ran::Domain;dir=Forward) branch(f,dom,ran,diff;dir=dir)


for TYP in (:FwdExpandingBranch,:RevExpandingBranch,:NeutralBranch)
  @eval @compat (b::$TYP)(x) = in(x,b.domain) ? unsafe_call(b,x) : throw(DomainError)
  @eval mapD(b::$TYP,x) = in(x,b.domain) ? unsafe_mapD(b,x) : error("DomainError: $x in $(b.domain)")
  @eval mapP(b::$TYP,x) = in(x,b.domain) ? unsafe_mapP(b,x) : throw(DomainError)
  @eval domain(b::$TYP) = b.domain
  @eval rangedomain(b::$TYP) = b.rangedomain
end
