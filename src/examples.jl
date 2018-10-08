export lanford, doubling, tupling

# Lanford map
for (fn,dfn,n) in ((:lanford_v0,:lanford_dv0,25),(:lanford_v1,:lanford_dv1,17))
  @eval ($fn)(x) = (5-sqrt($n-8x))/2
  @eval ($dfn)(x) = 2/sqrt($n-8x)
end
lanford_brk(T=Float64) = lanford_v0(one(T))

"""
    lanford(T=Float64)

Return the Lanford map, with type encoding T.

See also: [`tupling`](@ref), [`doubling`](@ref)
"""
lanford(T=Float64) = MarkovMap([lanford_v0,lanford_v1],[0..lanford_brk(T),lanford_brk(T)..1];dir=Reverse,diff=[lanford_dv0,lanford_dv1])

# tupling map
<<<<<<< HEAD
"""
    tupling(k::Int, d = Segment(0,1.))

Returns the full-branch interval map on domain d with k equally-sized branches.

See also: [`doubling`](@ref), [`lanford`](@ref)
"""
tupling(k::Int,d) = tupling(k,domain(d))
tupling(k::Int,d::Segment=Segment(0,1.)) = modulomap(x->k*(x-d.a)+d.a,d)
tupling(k::Int,d::PeriodicInterval) = FwdCircleMap(x->k*(x-d.a)+d.a,d,d,x->k*one(x))

"""
    doubling(d = Segment(0,1.))

Returns the full-branch interval map on domain d with 2 equally-sized branches.

See also: [`tupling`](@ref), [`lanford`](@ref)
"""
doubling(d=Segment(0,1.)) = tupling(2,d)
=======
tupling(k::Int,d=0..1.) = tupling(k,Domain(d))
tupling(k::Int,d::Segment) = modulomap(x->k*x,d,d,diff=x->oftype(x,k))#reversemodulomap(x->x/k,d,d,x->one(typeof(x))/k)
# tupling(k::Int,d::Segment) = MarkovMap([x->k*x-j for j = 0:k-1],[j/k..(j+1)/k for j = 0:k-1])
tupling(k::Int,d::PeriodicInterval) = RevCircleMap(x->x/k,d,d,x->one(typeof(x))/k)
doubling(d=0..1.) = tupling(2,d)
>>>>>>> 4c33d1554aa02b6a20399f0090f695de8ad3e517
