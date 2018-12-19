export rangedomain, Forward, Reverse

_NI(x) = error("not implemented: $x")

rangedomain(op::Operator) = domain(rangespace(op))

containstransfer(M::Operator) = ApproxFun.iswrapper(M::Operator) && containstransfer(M.op)
for optype in (:(ApproxFun.InterlaceOperator),:(ApproxFun.PlusOperator),:(ApproxFun.TimesOperator))
  eval(:(containstransfer(M::$optype) = any(containstransfer,M.ops)))
end

# Overloads

@compat const PureLaurent{D<:Domain} = Hardy{false,D}

@compat const GeneralInterval = Union{Segment,Interval,PeriodicSegment}

# rem()

# InterpolationNode

@compat struct InterpolationNode{S<:Space}<:Number
  sp::S
  k::Int
  n::Int
  function InterpolationNode{S}(sp::S,k::Int,n::Int) where S<:Space
    @assert k ≥ 1
    @assert k ≤ n
    new{typeof(sp)}(sp,k,n)
  end
end
InterpolationNode(sp::Space,k::Int,n::Int) = InterpolationNode{typeof(sp)}(sp,k,n)

Base.convert(::Type{T},p::InterpolationNode) where T<:Number = convert(T,points(p.sp,p.n)[p.k])
Base.convert(::Type{InterpolationNode{S}},p::InterpolationNode{S}) where S<:Space = p
Base.promote_rule(::Type{T},::Type{InterpolationNode{S}}) where {T<:Number,S<:Space} = promote_type(T,eltype(S))
Base.show(p::InterpolationNode) = show(convert(Nu))
ApproxFun.space(p::InterpolationNode) = p.sp
ApproxFun.domain(p::InterpolationNode) = domain(p.sp)
# eltype(p::InterpolationNode) = eltype(p.sp)
ApproxFun.prectype(p::InterpolationNode) = prectype(p.sp)

# Base.zero(p::InterpolationNode) = zero(eltype(p))

# Circle domains

interval_mod(x,a,b) = (b-a)*mod((x-a)/(b-a),1)+a
domain_mod(x,p::PeriodicSegment) = interval_mod(x,leftendpoint(p),rightendpoint(p))

# Interval domains

function approx_in(x,d::GeneralInterval; atol=50eps(arclength(d)), kwargs...)
    x <= rightendpoint(d)+atol && x >= leftendpoint(d) - atol
end
approx_in(x,d) = approx_in(convert(typeof(d),x),d)

function approx_issubset(a,b;kwargs...)
    acupb = a∪b
    atol_isapprox(acupb,b;kwargs...)
end

atol_isapprox(a,b; atol=50eps(max(arclength(a),arclength(b))), kwargs...) = isapprox(a,b;atol=atol,kwargs...)
coveringinterval(dsm::AbstractArray{T}) where {T<:AbstractInterval} = Interval(minimum(leftendpoint(convert(Domain,d)) for d in dsm),maximum(rightendpoint(convert(Domain,d)) for d in dsm))
coveringinterval(ds::AbstractArray) = coveringinterval([convert(Domain,d) for d in ds])
mapinterval(f,d::T) where T<:GeneralInterval = T(extrema((f(leftendpoint(d)),f(rightendpoint(d))))...)
mapinterval(f,d) = mapinterval(f,convert(Domain,d))

# Newton's method

domain_newton(f,df,y::U,D::Domain,x::U=convert(U,ApproxFun.checkpoints(D)[1]),tol=400eps(eltype(U))) where U = basic_newton(f,df,y,x,tol)
domain_guess(x,dom::Domain,ran::Domain) = rand(ran)

domain_guess(x::SVector{N},dom::ProductDomain,ran::ProductDomain) where N = SVector{N}([interval_guess(x[i],dom.domains[i],ran.domains[i]) for i =1:N])

function basic_newton(f,df,y::U,x::U,tol=40eps(U)) where U
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

# function monotonic_newton(f,df,y::U,xl::U,xh::U,tol=40eps(U)) where U ## TODO: REALLY WANT UPPER/LOWER GUESSES OR STH
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

function interval_newton(f,df,y::U,da::T,db::T,x::U=(da+db)/2,tol=40eps(max(abs(da),abs(db)))) where {T,U<:Real}
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

function interval_newton(f,df,y::U,da::T,db::T,x::U=(da+db)/2,tol=40eps(max(abs(da),abs(db)))) where {T,U<:Complex}
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
function domain_newton(f,df,y::U,D::GeneralInterval,x::U=(leftendpoint(D)+rightendpoint(D))/2,
            tol=40eps(max(leftendpoint(D),rightendpoint(D)))) where U<:Union{Real,Complex}
    # x = convert(U,x);
    interval_newton(f,df,y,leftendpoint(D),rightendpoint(D),x,tol)
end

interval_guess(y::Number,dom::Domain,ran::Domain) =
    (leftendpoint(dom)*rightendpoint(ran)-leftendpoint(ran)*rightendpoint(dom)+y*(rightendpoint(dom)-leftendpoint(dom)))/(rightendpoint(ran)-leftendpoint(ran))
domain_guess(y::Number,dom::GeneralInterval,ran::GeneralInterval) =
    interval_guess(y,dom,ran)

function disc_newton(f,df,y::U,rad::T,x::U=y,tol=40eps(rad)) where {T, U}
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

@compat struct FunctionDerivative{ff<:Function}
  f::ff
end
(fd::FunctionDerivative)(x) = oftype(x,dualpart(fd.f(Dual(x,one(x)))))

@compat struct BasisFun{S<:Space,II<:Integer}
  s::S
  k::II
end
(sk::BasisFun)(x) = getbasisfun(x,sk)
getbasisfun(x,sk::BasisFun) = Fun(sk.s,[zeros(prectype(sk.s),sk.k-1);one(prectype(sk.s))])(x)
getbasisfun_int(x,sk::BasisFun) = cumsum(Fun(sk.s,[zeros(prectype(sk.s),sk.k-1);one(prectype(sk.s))]))(x)

@compat struct BasisFunInt{S<:Space,II<:Integer}
  sk::BasisFun{S,II}
end
cumsum(sk::BasisFun) = BasisFunInt(sk)
(ski::BasisFunInt)(x) = getbasisfun_int(ski.sk,x)

#specialised getbasisfuns

fourierCSk(x,d,k::Integer) = rem(k,2) == 1 ? cos(fld(k,2)*tocanonical(d,x)) : sin(fld(k,2)*tocanonical(d,x))
function fourierCSk_int(x,d,k::Integer)
  k == 1 && return x
  (rem(k,2) == 1 ? -sin(fld(k,2)*tocanonical(d,x)) : cos(fld(k,2)*tocanonical(d,x)))/
    fld(k,2)/tocanonicalD(d,x)
end
getbasisfun(x,sk::BasisFun{Fourier{DD},K}) where {DD, K<:Integer} = fourierCSk(x,domain(sk.s),sk.k)
getbasisfun_int(x,sk::BasisFun{Fourier{DD},K}) where {DD, K<:Integer} = fourierCSk_int(x,domain(sk.s),sk.k)

chebyTk(x,d,k::Integer) = cos((k-1)*acos(tocanonical(d,x))) #roundoff error grows linearly(??) with k may not be bad wrt x too
function chebyTk_int(x,d,k::Integer)
#   k == 1 && return tocanonical(d,x)/tocanonicalD(d,x)
  k == 2 && return tocanonical(d,x)^2/2tocanonicalD(d,x)
  ((k-1)*chebyTk(x,d,k+1) - k*tocanonical(d,x)*chebyTk(x,d,k))/((k-1)^2-1)/tocanonicalD(d,x)
end
getbasisfun(x,sk::BasisFun{F,K}) where {F<:Chebyshev, K<:Integer} = chebyTk(x,domain(sk.s),sk.k)
getbasisfun_int(x,sk::BasisFun{F,K}) where {F<:Chebyshev,K<:Integer} = chebyTk_int(x,domain(sk.s),sk.k)

function getbasisfun(x,sk::BasisFun{F,K}) where {F<:TensorSpace,K<:Integer}
  ks = ApproxFun.tensorizer(sk.s)[sk.k]
  prod(getbasisfun(x[i],BasisFun(sk.s.spaces[i],ks[i])) for i = eachindex(sk.s.spaces))
end
# no getbasisfun_int as you don't have antiderivatives


# # Directions
@compat abstract type Direction; end
@compat struct ForwardDirection <: Direction; end
@compat struct ReverseDirection <: Direction; end
const Forward = ForwardDirection()
const Reverse = ReverseDirection()
