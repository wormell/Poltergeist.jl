function transferfunction{M,B}(x,m::InducedMarkovMap{M,B},d::Space,k::Integer,T)
  v = convert(T,x)
  dv = one(T)
  y = zero(T)

  while abs(dv) > eps(T)
    y+= abs(dv)*transferfunction(v,m.m,d,k,T)
    a = mapinvP(m.b,v)
    v = a[1]
    dv *= a[2]
  end

  y
end

function transferfunction{M,B<:NeutralBranch}(x,m::InducedMarkovMap{M,B},d::Space,k::Integer,T)
  v = convert(T,x)
  dv = one(T)
  y = zero(T)

  while abs(toring(m.b.ab,v)) < m.b.ab.rd
  #for i = 1:K
    y+= abs(dv)*transferfunction(v,m.m,d,k,T)
    a = mapinvP(m.b,v)
    v = a[1]
    dv *= a[2]
  end

  y += dv*fullmeasure_asym(m.b.ab,m.m,v,d,k,T)

  y
end




function fullmeasure_asym(R::AbelFunction,m::MarkovMap,x::Number,d::Space,k::Integer,T)
#  @assert abs(toring(R,x)) > R.rd
  IQ = real(mapasymD(R,x)) * transferfunction_int(R.p,x,m,d,k,T)
  JQ = transferfunction(x,m,d,k,T)/2

  #CQ
  K = convert(Int,ceil(log2(R.alpha * R.dh0 / 2pi * -log(eps(T)))))
  N = 2^convert(Int,ceil(log2(-log(eps(T))/log(1+sqrt(2)))))
#  edges = [0.;2.0.^(0:K)]
  mids = T[1;3*(2).^(0:K-1)]/2
  scales = T[1;(2).^(0:K-1)]/2
  intvals = zeros(T,N)
  mα = -R.alpha
  transfer_x = transferfunction(x,m,d,k,T)
  t2 = 90(eps(T))^(1//3)*-mα*x^mα
  for (i,s) in enumerate(chebyshevpoints(T,N))
    for j in eachindex(mids)
      t = mids[j]+scales[j]*s
      ω = (x^mα + im*t)^(1/mα)
#       Δ = -pi*(abs(t*x^(-mα)/mα) < 20sqrt(eps(T)) ? -im*t*mapasymD(R,x)/R.alpha*x^(1-mα) + t^2*x^(1-2mα)/2mα^2 *(-(1-mα)/mapasymD(R,x) + x*dualpart(mapasymD(R,Dual(x,one(x))))):
#               (mapasym(R,x)-mapasym(R,ω)))
      if t > t2
        Δ = -pi*(mapasym(R,x)-mapasym(R,ω))
      else
        ω2 = (x^mα + im*t2)^(1/mα)
        Δ = -pi*(mapasym(R,x)-mapasym(R,ω2))
        Δ = real(Δ)*(t/t2)^2 + im *imag(Δ)*(t/t2)
      end
      τ = abs(Δ) < sqrt(eps(T)) ? cot(Δ)-im : 2im/(exp(2im*Δ)-1)
#       i == N && println((t,abs(τ),abs(intvals[i])))
      intvals[i] += scales[j]*real((transferfunction(ω,m,d,k,Complex{T})*ω^(1-mα)-transfer_x*x^(1-mα)) * τ)
      # putting back in the transfer_x*x^(1-mα)) term
      tanx = tan(real(Δ));
      intvals[i] += scales[j]*transfer_x *x^(1-mα) * tanx/(sinh(imag(Δ))^2+tanx^2*cosh(imag(Δ))^2)
    end
    #     println((t,abs(w-x),dR,intvals[i]))
  end
 #   println(intvals)
  dL = ApproxFun.Chebyshev()
  CQ = real(mapasymD(R,x))/R.alpha * sum(Fun(dL,ApproxFun.transform(dL,intvals)))

  IQ+JQ+CQ
end

type FullAcim{A,B,T,D,R,MM<:InducedMarkovMap}
  L::CachedOperator{A,B,ConcreteTransfer{T,D,R,MM}}
  rho::Fun{T,D}
end
FullAcim{A,B,T,D,R,MM<:InducedMarkovMap}(L::CachedOperator{A,B,ConcreteTransfer{T,D,R,MM}}) = FullAcim(L,acim(L))

@compat function (rhof::FullAcim)(x)
  dot(rhof.rho.coefficients,eltype(rhof.L)[transferfunction(x,getmap(rhof.L),domainspace(rhof.L.op),i,eltype(L)) for i = 1:length(rhof.rho.coefficients)])
end
