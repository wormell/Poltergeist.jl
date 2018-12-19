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
"""
    tupling(k::Int, d = Interval(0,1.))

Returns the full-branch interval map on domain d with k equally-sized branches.

See also: [`doubling`](@ref), [`lanford`](@ref)
"""
tupling(k::Int,d=Interval(0,1.)) = modulomap(x->k*(x-leftendpoint(d))+leftendpoint(d),d)
tupling(k::Int,d::PeriodicSegment) = RevCircleMap(x->(x-leftendpoint(d))/k+leftendpoint(d),d,d,x->one(x)/k)

"""
    doubling(d = Interval(0,1.))

Returns the full-branch interval map on domain d with 2 equally-sized branches.

See also: [`tupling`](@ref), [`lanford`](@ref)
"""
doubling(d=Interval(0,1.)) = tupling(2,d)
