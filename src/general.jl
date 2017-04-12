export rangedomain
rangedomain(op::Operator) = domain(rangespace(op))

containstransfer(M::Operator) = ApproxFun.iswrapper(M::Operator) && containstransfer(M.op)
for optype in (:(ApproxFun.InterlaceOperator),:(ApproxFun.PlusOperator),:(ApproxFun.TimesOperator))
  eval(:(containstransfer(M::$optype) = any(containstransfer,M.ops)))
end

# Overloads

typealias PureLaurent{D<:Domain} Hardy{false,D}

typealias GeneralInterval Union{Segment,PeriodicInterval}

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

# Circle domains

interval_mod(x,a,b) = (b-a)*mod((x-a)/(b-a),1)+a
mod(x,p::PeriodicInterval) = interval_mod(x,p.a,p.b)

# Newton's method

domain_newton{U}(f,df,y::U,D::Domain,x::U=convert(U,rand(D)),tol=400eps(eltype(U))) = basic_newton(f,df,y,x,tol)
domain_guess(x,dom::Domain,ran::Domain) = rand(ran)

domain_guess{N}(x::Vec{N},dom::ProductDomain,ran::ProductDomain) = Vec{N}([interval_guess(x[i],dom.domains[i],ran.domains[i]) for i =1:N])

function basic_newton{U}(f,df,y::U,x::U,tol=40eps(U))
  z = x
  rem = f(z)-y
  for i = 1:2log2(-200log(eps(eltype(y))))
    #    while abs(rem) > 30eps(maximum(∂(D)))
    z -= df(z) \ rem
    rem = f(z)-y
    norm(rem) < tol && break
  end

  norm(rem) > tol && error("Newton: failure to converge: y = $y, x estimate = $z, rem = $rem")
  z
end

# function monotonic_newton{U}(f,df,y::U,xl::U,xh::U,tol=40eps(U)) ## TODO: REALLY WANT UPPER/LOWER GUESSES OR STH
#   σ = sign(df(x))
#   det(σ) == 0 && error("f has critical point at $x")
#   z = xl
#   σzl = σ*xl; σzu = σ*xh; # actually find x,u
#   rem = f(z)-y
#   σz = σ*z
#   rem < 0 ? (σzu = min(σzu,σz)) : (σzl = max(σzl,σz))
#   println((z,y,rem,σz,σzl,σzu))
#   for i = 1:200#log2(-200log(eps(eltype(y))))
#     #    while abs(rem) > 30eps(maximum(∂(D)))
#     z -= df(z) \ rem
#     rem = f(z)-y
#     norm(rem) < tol && break
#
#     println((z,y,rem,σz,σzl,σzu))
#
#     σz = σ*z
#     rem < 0 ? (σzu = min(σzu,σz)) : (σzl = max(σzl,σz))
#     σzu < σzl && (σzu += 4tol, σzl -= 4tol)
#     (σz > σzu || σz < σzl) && (z = σ*(σzl + (σzu-σzl)*rand()))
#     println((z,y,rem,σz,σzl,σzu))
#   end
#
#   norm(rem) > tol && error("Newton: failure to converge")
#   z
# end

function interval_newton{T,U<:Real}(f,df,y::U,da::T,db::T,x::U=(da+db)/2,tol=40eps(max(abs(da),abs(db))))
  rem = f(x)-y
  for i = 1:2log2(-200log(eps(typeof(y))))
    #    while abs(rem) > 30eps(maximum(∂(D)))
    x -= rem / df(x) + eps(x)*(mod(i,11)-5)/4
    x .> db && (x = db)
    x .< da && (x = da)
    rem = f(x)-y
    #  println((x,f(x)-y))
    abs(rem) < tol && break
  end

  abs(rem) > tol && error("Newton: failure to converge: y = $y, x estimate = $x, rem = $rem")
  x
end

function interval_newton{T,U<:Complex}(f,df,y::U,da::T,db::T,x::U=(da+db)/2,tol=40eps(max(abs(da),abs(db))))
  rem = f(x)-y
  for i = 1:2log2(-200log(eps(typeof(abs(y)))))
    #    while abs(rem) > 30eps(maximum(∂(D)))
    x -= rem / df(x)
    rem = f(x)-y
    abs(rem) < tol && break
  end

  abs(rem) > tol && error("Newton: failure to converge: y = $y, x estimate = $x, rem = $rem")
  x
end
domain_newton{U<:Union{Real,Complex}}(f,df,y::U,D::GeneralInterval,x::U=(D.a+D.b)/2,tol=40eps(max(abs(D.a),abs(D.b)))) = interval_newton(f,df,y,D.a,D.b,x,tol)

interval_guess(y::Number,dom::Domain,ran::Domain) = (dom.a*ran.b-ran.a*dom.b+y*(dom.b-dom.a))/(ran.b-ran.a)
domain_guess(y::Number,dom::GeneralInterval,ran::GeneralInterval) = interval_guess(y,dom,ran)

function disc_newton{T,U}(f,df,y::U,rad::T,x::U=y,tol=40eps(rad))
#   x = convert(typeof(y),rad*rand(typeof(real(y))))
  rem = f(x) - y
  for i = 1:2log2(-200log(eps(typeof(real(y)))))
    x -= rem / df(x)
    rem = f(x) - y
    abs(rem) < tol && break
    abs(x) > rad && (x /= abs(x))
  end
  abs(rem) > tol && error("Newton: failure to converge: y = $y, x estimate = $x, rem = $rem")
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
@compat (fd::FunctionDerivative)(x) = oftype(x,dualpart(fd.f(Dual(x,one(x)))))

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
getbasisfun{DD,K<:Integer}(x,sk::BasisFun{Fourier{DD},K},T) = fourierCSk(x,domain(sk.s),sk.k)
getbasisfun_int{DD,K<:Integer}(x,sk::BasisFun{Fourier{DD},K},T) = fourierCSk_int(x,domain(sk.s),sk.k)

chebyTk(x,d,k::Integer) = cos((k-1)*acos(tocanonical(d,x))) #roundoff error grows linearly(??) with k may not be bad wrt x too
function chebyTk_int(x,d,k::Integer)
#   k == 1 && return tocanonical(d,x)/tocanonicalD(d,x)
  k == 2 && return tocanonical(d,x)^2/2tocanonicalD(d,x)
  ((k-1)*chebyTk(x,d,k+1) - k*tocanonical(d,x)*chebyTk(x,d,k))/((k-1)^2-1)/tocanonicalD(d,x)
end
getbasisfun{F<:Chebyshev,K<:Integer}(x,sk::BasisFun{F,K},T) = chebyTk(x,domain(sk.s),sk.k)
getbasisfun_int{F<:Chebyshev,K<:Integer}(x,sk::BasisFun{F,K},T) = chebyTk_int(x,domain(sk.s),sk.k)

function getbasisfun{F<:TensorSpace,K<:Integer}(x,sk::BasisFun{F,K},T)
  ks = ApproxFun.tensorizer(sk.s)[sk.k]
  @compat prod(getbasisfun(x[i],BasisFun(sk.s.spaces[i],ks[i]),T) for i = eachindex(sk.s.spaces))
end
# no getbasisfun_int as you don't have antiderivatives


# # Directions
# abstract FunctionDirection
# type Fwd <: FunctionDirection; end
# type Rev <: FunctionDirection; end
