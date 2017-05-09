function transferfunction{M,B}(x,m::InducedMarkovMap{M,B},sk,T)
  v = convert(T,x)
  dv = one(T)
  y = zero(T)

  while abs(dv) > eps(T)
    y+= abs(dv)*transferfunction(v,m.m,sk,T)
    a = mapinvP(m.b,v)
    v = a[1]
    dv *= a[2]
  end

  y
end
