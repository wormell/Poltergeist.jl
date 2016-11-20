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
space(p::InterpolationNode) = p.sp
domain(p::InterpolationNode) = domain(p.sp)
Base.eltype(p::InterpolationNode) = eltype(p.sp)

Base.zero(p::InterpolationNode) = zero(eltype(p))

# Newton's method

function interval_newton{T}(f,df,y::Number,da::T,db::T,tol=10eps(max(abs(da),abs(db))))
  x = da+(db-da)*rand(typeof(y))
  rem = f(x)-y
  for i = 1:2log2(-200log(eps(typeof(y))))
    #    while abs(rem) > 30eps(maximum(∂(D)))
    x -= rem / df(x)
    rem = f(x)-y
    abs(rem) < tol && break
    x .> db && (x = db)
    x .< da && (x = da)
  end
  abs(rem) > tol && error("Newton: failure to converge")
  x
end
interval_newton(f,df,y::Number,D::Domain,tol=10eps(max(abs(D.a),abs(D.b)))) = interval_newton(f,df,y,D.a,D.b,tol)

function disc_newton{T}(f,df,y::Number,rad::T,tol=10eps(rad))
  x = rad*rand()
  rem = f(x) - y
  for i = 1:2log2(-200log(eps(typeof(y))))
    x -= rem / df(x)
    rem = f(x) - y
    abs(rem) < tol && break
    abs(x) > rad && (x /= abs(x))
  end
  abs(rem) > tol && error("Newton: failure to converge")
  x
end

# ForwardDiff
function forwarddiff(f)
  function dfdx(x)
    ForwardDiff.derivative(f,x)::typeof(x)
  end
  dfdx
end
