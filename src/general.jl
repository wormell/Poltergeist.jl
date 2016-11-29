export rangedomain
rangedomain(op::Operator) = domain(rangespace(op))

containstransfer(M::Operator) = ApproxFun.iswrapper(M::Operator) && containstransfer(M.op)
for optype in (:(ApproxFun.InterlaceOperator),:(ApproxFun.PlusOperator),:(ApproxFun.TimesOperator))
  eval(:(containstransfer(M::$optype) = any(containstransfer,M.ops)))
end

# Overloads

typealias PureLaurent{D<:Domain} Hardy{false,D}

#ApproxFun.qrfact(A::ApproxFun.QROperator) = A
# rem()

# InterpolationNode

type InterpolationNode{S<:Space}<:Number
  sp::S
  k::Int
  n::Int
  function InterpolationNode(sp::Space,k::Int,n::Int)
    @assert k ≥ 1
    @assert k ≤ n
    new{typeof(sp)}(sp,k,n)
  end
end
InterpolationNode(sp::Space,k::Int,n::Int) = InterpolationNode{typeof(sp)}(sp,k,n)

Base.convert{T<:Number}(::Type{T},p::InterpolationNode) = convert(T,points(p.sp,p.n)[p.k])
Base.convert{S<:Space}(::Type{InterpolationNode{S}},p::InterpolationNode{S}) = p
Base.promote_rule{T<:Number,S<:Space}(::Type{T},::Type{InterpolationNode{S}}) = promote_type(T,eltype(S))
Base.show(p::InterpolationNode) = show(convert(Nu))
ApproxFun.space(p::InterpolationNode) = p.sp
ApproxFun.domain(p::InterpolationNode) = domain(p.sp)
Base.eltype(p::InterpolationNode) = eltype(p.sp)

Base.zero(p::InterpolationNode) = zero(eltype(p))

# Newton's method

function interval_newton{T,U<:Real}(f,df,y::U,da::T,db::T,x::U=(da+(db-da)*rand(typeof(y))),tol=10eps(max(abs(da),abs(db))))
  rem = f(x)-y
  for i = 1:2log2(-200log(eps(typeof(y))))
    #    while abs(rem) > 30eps(maximum(∂(D)))
    x -= rem / df(x)
    x .> db && (x = db)
    x .< da && (x = da)
    rem = f(x)-y
    abs(rem) < tol && break
  end

  abs(rem) > tol && error("Newton: failure to converge")
  x
end
interval_newton(f,df,y::Number,D::Domain,tol=10eps(max(abs(D.a),abs(D.b)))) = interval_newton(f,df,y,D.a,D.b,tol)
interval_guess(y::Number,dom::Domain,ran::Domain) = (dom.a*ran.b-ran.a*dom.b+y*(dom.b-dom.a))/(ran.b-ran.a)

function interval_newton{T,U<:Complex}(f,df,y::U,da::T,db::T,x::U=(da+(db-da)*rand(typeof(y))),tol=10eps(max(abs(da),abs(db))))
  rem = f(x)-y
  for i = 1:2log2(-200log(eps(typeof(abs(y)))))
    #    while abs(rem) > 30eps(maximum(∂(D)))
    x -= rem / df(x)
    rem = f(x)-y
    abs(rem) < tol && break
  end

  abs(rem) > tol && error("Newton: failure to converge")
  x
end

function disc_newton{T,U}(f,df,y::U,rad::T,x::U=y,tol=10eps(rad))
#   x = convert(typeof(y),rad*rand(typeof(real(y))))
  rem = f(x) - y
  for i = 1:2log2(-200log(eps(typeof(real(y)))))
    x -= rem / df(x)
    rem = f(x) - y
    abs(rem) < tol && break
    abs(x) > rad && (x /= abs(x))
  end
  abs(rem) > tol && error("Newton: failure to converge")
  x
end

# ForwardDiff
# function forwarddiff(f)
#   function dfdx(x)
#     ForwardDiff.derivative(f,x)::typeof(x)
#   end
#   dfdx
# end

type FunctionDerivative{ff<:Function}
  f::ff
end
@compat (fd::FunctionDerivative)(x::Number) = oftype(x,dualpart(fd.f(Dual(x,one(x)))))

type BasisFun{S<:Space,II<:Integer}
  s::S
  k::II
end
@compat (b::BasisFun)(x) = getbasisfun(x,sk,eltype(sk.s))
getbasisfun(x,sk::BasisFun,T) = Fun(sk.s,[zeros(T,sk.k-1);one(T)])(x)
getbasisfun_int(x,sk::BasisFun,T) = cumsum(Fun(sk.s,[zeros(T,sk.k-1);one(T)]))(x)

#specialised getbasisfuns

fourierCSk(x,d,k::Integer) = rem(k,2) == 1 ? cos(fld(k,2)*tocanonical(d,x)) : sin(fld(k,2)*tocanonical(d,x))
function fourierCSk_int(x,d,k::Integer)
  k == 1 && return x
  (rem(k,2) == 1 ? -sin(fld(k,2)*tocanonical(d,x)) : cos(fld(k,2)*tocanonical(d,x)))/
    fld(k,2)/tocanonicalD(d,x)
end
getbasisfun{F<:Fourier,K<:Integer}(x,sk::BasisFun{F,K},T) = fourierCSk(x,domain(sk.s),sk.k)
getbasisfun_int{F<:Fourier,K<:Integer}(x,sk::BasisFun{F,K},T) = fourierCSk_int(x,domain(sk.s),sk.k)

chebyTk(x,d,k::Integer) = cos((k-1)*acos(tocanonical(d,x))) #roundoff error grows linearly(??) with k may not be bad wrt x too
function chebyTk_int(x,d,k::Integer)
#   k == 1 && return tocanonical(d,x)/tocanonicalD(d,x)
  k == 2 && return tocanonical(d,x)^2/2tocanonicalD(d,x)
  ((k-1)*chebyTk(x,d,k+1) - k*tocanonical(d,x)*chebyTk(x,d,k))/((k-1)^2-1)/tocanonicalD(d,x)
end
getbasisfun{F<:Chebyshev,K<:Integer}(x,sk::BasisFun{F,K},T) = chebyTk(x,domain(sk.s),sk.k)
getbasisfun_int{F<:Chebyshev,K<:Integer}(x,sk::BasisFun{F,K},T) = chebyTk_int(x,domain(sk.s),sk.k)
