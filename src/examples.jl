export lanford, doubling, tupling

# Lanford map
for (fn,dfn,n) in ((:lanford_v0,:lanford_dv0,25),(:lanford_v1,:lanford_dv1,17))
  @eval ($fn)(x) = (5-sqrt($n-8x))/2
  @eval ($dfn)(x) = 2/sqrt($n-8x)
end
lanford_brk(T=Float64) = lanford_v0(one(T))
lanford(T=Float64) = MarkovMap([lanford_v0,lanford_v1],[0..lanford_brk(T),lanford_brk(T)..1];dir=Reverse,diff=[lanford_dv0,lanford_dv1])

# tupling map
tupling(k::Int,d=0..1.) = tupling(k,Domain(d))
tupling(k::Int,d::Segment) = modulomap(x->k*x,d,d,diff=x->oftype(x,k))#reversemodulomap(x->x/k,d,d,x->one(typeof(x))/k)
tupling(k::Int,d::PeriodicInterval) = RevCircleMap(x->x/k,d,d,x->one(typeof(x))/k)
doubling(d) = tupling(2,d)
