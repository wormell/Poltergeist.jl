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

#   function AbelFunction(h,dh,alpha,p,sgn,r0)
#     new(h,dh,Complex{T}[],alpha,p,r0,sgn,zero(T),zero(T),dh(0))
#   end
end
AbelFunction(h,dh,cfs,alpha,p,sgn,r0,rd,rab,dh0) = AbelFunction{typeof{h},typeof{dh}}

toring(R::AbelFunction,x) = -(R.sgn*(x-R.p)).^(-R.alpha)/(R.alpha*R.dh0)
toringD(R::AbelFunction,x) = -(R.sgn*(x-R.p)).^(-R.alpha)/(x-R.p)
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


function construct_asymptotic_approx!{T<:Real}(R::AbelFunction{T},r0::T=R.r0)
  N = convert(Int,ceil(-log(eps(T))/2))+5
#  R.cfs = zeros(T,N)

  Qr0 = Q(R)
  Cs_initialguess = Cs(R)
  R.rd = max(R.rd,N*Qr0*(1+Qr0*exp(-1/N)),(2Cs_initialguess+1)*ringhconv(R,r0))
#  println(R.rd)
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



function AbelFunction{T<:Real,ff,gg}(h::ff,dh::gg,r0::T,alpha::T,p::T=zero(T),sgn::T=one(T))
  R=AbelFunction(h,dh,Complex{T}[],alpha,p,sgn,r0,zero(T),zero(T),convert(T,dh(0)))
  construct_asymptotic_approx!(R)
  R
end

function mapasym(R::AbelFunction,x)
  #  if abs(x) < fromring(R.rd)
  y = zero(x)
  z = toring(R,x)
  for i = length(R.cfs):-1:2
    y += R.cfs[i]*z.^(1-i)/(1-i)
  end
  y += R.cfs[1]*log(-z)
  y += z
  y
#   else
#     error("Absolute value of x too big")
#   end
end

function (R::AbelFunction)(x::Number)
  if abs(toring(R,x)) > R.rd
    return mapasym(R,x)
  else
    n = 0
    y = ringhconv(R,toring(R,x))
    println((x,y,ringhconv(R,R.rd)))
    while abs(y) > abs(ringhconv(R,R.rd))
      y = disc_newton(gg->(gg*(R.h(gg))^R.alpha)::typeof(gg),gg->((gg*R.alpha*R.dh(gg)+R.h(gg))*R.h(gg)^(R.alpha-1))::typeof(gg),y,max(2y,y/(1-y)))
      n += 1
      n == 100 && error("Backwards iteration failed in Abel function")
    end
    println(n)
    return mapasym(R,fromring(R,ringhconv(R,y))) + n
  end
end
