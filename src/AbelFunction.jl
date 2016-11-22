export AbelFunction

type AbelFunction{T<:Real,ff,gg}
  h::ff #h
  dh::gg #derivative of h
  cfs::Vector{Complex{T}} #coeffs of asymptotic approximation to A'
  alpha::T # power
  p::T #fixed point
  sgn::T #which way neutral fixed point goes
  r0::T #radius of analyticity of h
  rd::T #radius of A
  rab::T #radius of A inverse??
  #  K::Int #
  dh0::T # first derivative of h at 0
  offset::T # (real) constant offset of coefficients

  function AbelFunction(h,dh,cfs,alpha,p,sgn,r0,rd,rab,dh0,offset)
    @assert alpha > 0
    @assert sgn^2 == one(T)
    @assert r0 > 0
    @assert rd ≥ 0
    @assert rab ≥ 0
    @assert dh0 ≥ 0
    new(h,dh,cfs,alpha,p,sgn,r0,rd,rab,dh0,offset)
  end
end
function AbelFunction{T<:Real}(h,dh,cfs::Vector{Complex{T}},alpha::T,p::T,sgn::T,r0::T,rd::T,rab::T,dh0::T,offset::T)
  AbelFunction{T,typeof(h),typeof(dh)}(h,dh,cfs,alpha,p,sgn,r0,rd,rab,dh0,offset)
end

toring(R::AbelFunction,x) = -(R.sgn*(x-R.p)).^(-R.alpha)/(R.alpha*R.dh0)
toringD(R::AbelFunction,x) = (R.sgn*(x-R.p)).^(-R.alpha)/(x-R.p)/R.dh0
fromring(R::AbelFunction,z) = R.p + R.sgn*(-R.alpha*R.dh0*z).^(-1/R.alpha)
fromringD(R::AbelFunction,z) = R.sgn*R.dh0*(-R.alpha*R.dh0*z).^(-1/R.alpha - 1)

ringhconv(R::AbelFunction,y) = -1/(R.alpha*R.dh0*y)
ringhconvD(R::AbelFunction,y) = 1/(R.alpha*R.dh0*y^2)
#hdomtoring(R::AbelFunction,z) = -1/()


function Q(R::AbelFunction,rlaur::Real=abs(ringhconv(R,R.r0)))
  hd{T}(z::T) = (R.h(ringhconv(R,z)).^R.alpha - 1)::T
  rlaur*sum(abs(real(coefficients(Fun(hd,PureLaurent(Circle(rlaur)))))))
end
function Cs(R::AbelFunction,s::Real=abs(ringhconv(R,R.r0)))
  hd{T}(z::T) = (R.h(ringhconv(R,z)).^(-R.alpha) - 1 - 1/z)::T
  sum(abs(coefficients(Fun(hd,PureLaurent(Circle(s))))))
end


function constructasym!{T<:Real}(R::AbelFunction{T},r0::T=R.r0)
  N = convert(Int,ceil(-log(eps(T))/2))+5
  #  R.cfs = zeros(T,N)

  Qr0 = Q(R)
  Cs_initialguess = Cs(R)
  R.rd = max(R.rd,N*Qr0*(1+Qr0*exp(-1/N)),(2Cs_initialguess+1)*ringhconv(R,r0))
  rd_space = PureLaurent(Circle(R.rd))

  Delta = zeros(T,N,N)

  funcfs = coefficients(Fun(z->-R.alpha*log(R.h(ringhconv(R,z)))::typeof(z),rd_space)')
  Delta[:,1] = real(pad!(funcfs[2:end],N))
  for k = 2:N
    funcfs = coefficients(Fun(z->((R.h(ringhconv(R,z))^(R.alpha*(k-1))-1)/(1-k)*z^(1-k))::typeof(z),rd_space)')
    Delta[k:end,k] = real(pad!(funcfs[k+1:end],N-k+1))
  end
  #Delta[abs(Delta) .< 10*N*eps(T)] = zero(T)
  #  R.cfs = Delta \ real(pad!(coefficients(Fun(z->(z/R.h(ringhconv(R,z))^(R.alpha)-z-1)::typeof(z),rd_space)')[2:end],N))
  acfs=-coefficients(Fun(z->(z/R.h(ringhconv(R,z))^(R.alpha)-z-1)::typeof(z),rd_space)')
  Rcfs = (Delta\real(pad!(acfs[2:end],N)))./R.rd.^(1:N)
  chop!(Rcfs,N*eps(T))
  R.cfs = Rcfs .* R.rd.^(1:length(Rcfs))

  R
end

function AbelFunction{T<:Real,ff,gg}(h::ff,dh::gg,alpha::T,r0::T,p::T=zero(T),sgn::T=one(T))
  R=AbelFunction(h,dh,Complex{T}[],alpha,p,sgn,r0,zero(T),zero(T),convert(T,dh(0)),zero(T))
  constructasym!(R)
  R
end

hdisc_guess(R::AbelFunction,y) = abs(ringhconv(R,y)) > 3 ? ringhconv(R,ringhconv(R,y)-1) : y

function hdisc_newton{T,U}(R::AbelFunction{T},y::U,rad::T=abs(y)/max(0.1,1-abs(y)),
                         x::U=hdisc_guess(R,y),tol=10eps(rad))
  x = convert(typeof(y),rad*rand(typeof(real(y))))
  rhx = convert(T,R.h(x))
  rem = x*rhx^R.alpha - y
  for i = 1:2log2(-200log(eps(one(real(y)))))
    x -= rem / rhx^R.alpha / (1 + R.alpha*x*R.dh(x)/rhx)
    rhx = R.h(x)
    rem = x*(R.h(x))^R.alpha  - y
    abs(rem) < tol && break
    abs(x) > rad && (x /= abs(x))
  end
  abs(rem) > tol && error("Newton: failure to converge")
  x

end

function mapasym(R::AbelFunction,x)
  y = zero(x)
  z = toring(R,x)
  for i = length(R.cfs):-1:2
    y += R.cfs[i]*z.^(1-i)/(1-i)
  end
  y += R.cfs[1]*log(-z)
  y += R.offset
  y += z
  y
end
function mapasymD(R::AbelFunction,x)
  y = zero(x)
  z = toring(R,x)
  for i = length(R.cfs):-1:1
    y += R.cfs[i]*z.^(-i)
  end
  y += 1
  y*toringD(R,x)
end
mapasymP(R::AbelFunction,x) = (mapasym(R,x),mapasymD(R,x))

function mapP{T}(R::AbelFunction{T},x::Number)
  if abs(toring(R,x)) > R.rd
    return mapasymP(R,x)
  else
    n = 0
    y = ringhconv(R,toring(R,x))
    dv = 1/(toringD(R,x)*ringhconvD(R,toring(R,x)))
    while abs(y) > abs(ringhconv(R,R.rd))
      y = hdisc_newton(R,y)
      rhy = R.h(y)
      dv *= rhy^R.alpha * (1 + R.alpha*y*R.dh(y)/rhy)
      n += 1
      n == 200 && error("Backwards iteration failed in Abel function")
    end
    return (mapasym(R,fromring(R,ringhconv(R,y))) + n,mapasymD(R,fromring(R,ringhconv(R,y)))/dv)
  end
end
(R::AbelFunction)(x::Number) = mapP(R,x)[1]
mapD(R::AbelFunction,x::Number) = mapP(R,x)[2]


function zeroat!{T}(R::AbelFunction{T},x::Number)
  R.offset -= convert(T,R(x))
  R
end

# TODO: Abel function inverse
#mapasym_recip(R::AbelFunction,x::Number) =1/mapasym(R,x)
#mapasym_recipD(R::AbelFunction,x::Number) = mapasymD(R,x)/mapasym(R,x)^2
#function mapasyminv(R::AbelFunction,x::Number)


