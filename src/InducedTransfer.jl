function transferfunction(x,m::InducedMarkovMap,d::Space,k::Integer,T)
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

function transferfunction{D,R,M,B<:NeutralBranch}(x,m::InducedMarkovMap{D,R,M,B},d::Space,k::Integer,T)
  v = convert(T,x)
  dv = one(T)
  y = zero(T)

  while abs(toring(m.b.ab.r,v)) > m.b.ab.rd
    y+= abs(dv)*transferfunction(x,m.m,d,k,T)
    a = mapinvP(L.m.b,v)
    v = a[1]
    dv *= a[2]
  end

  y += dv*fullmeasure_asym(m.b.ab,m.m,v,d,k,T)

  y
end




function fullmeasure_asym(R::AbelFunction,M::MarkovMap,x::Number,d::Space,k::Integer,T)
  IQ = mapD(R,x) * transferfunction_int(R.p,x,m,d,k,T)
  JQ = transferfunction(x,m,d,k,T)/2
  #CQ ... heh
  IQ + JQ
end