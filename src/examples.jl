export lanford, doubling, tupling

# Lanford map
for (fn,dfn,n) in ((:lanford_v0,:lanford_dv0,25),(:lanford_v1,:lanford_dv1,17))
  @eval ($fn)(x) = (5-sqrt($n-8x))/2
  @eval ($dfn)(x) = 2/sqrt($n-8x)
end
lanford_brk(T=Float64) = lanford_v0(one(T))
lanford(T=Float64) = MarkovMap([lanford_v0,lanford_v1],[0..lanford_brk(T),lanford_brk(T)..1];dir=Reverse,diff=[lanford_dv0,lanford_dv1])

# tupling map
tupling(k::Int,d) = tupling(k,domain(d))
tupling(k::Int,d::Segment=0..1) = modulomap(x->k*x,d)
tupling(k::Int,d::PeriodicInterval) = FwdCircleMap(x->k*x,d,d,x->k*one(x))
doubling(d) = tupling(2,d)
