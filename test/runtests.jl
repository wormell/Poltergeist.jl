using GhostCoop
using Base.Test
using ApproxFun

# write your own tests here

# Expanding tests

f1(x)=2x+sin(2pi*x)/8pi; f2(x)=2x+sin(2pi*x)/8pi-1
fv1(x) = x/2+sin(2pi*x)/8pi; fv2(x) = x/2+1/2+sin(2pi*x)/8pi
fv1d(x) = 1/2+cos(2pi*x)/4; fv2d(x) = 1/2+cos(2pi*x)/4

# Periodic domain
d1 = PeriodicInterval([0,1])
M1 = (GC.MarkovMap(d1,fv,[0,0.5,1],"rev"));
GC.acim(M1)
@time rho1 = GC.acim(M1)

M1a = GC.MarkovMap(d1,[f1,f2],[f1d,f2d],[0,0.5,1])
GC.acim(M1a)
@time rho1a = GC.acim(M1a)

# Non-periodic domain
d2 = Interval([0,1])
M2 = GC.MarkovMap(d2,fv,[0,0.5,1],"rev");
GC.acim(M2)
@time rho2 = GC.acim(S2)

M2a = GC.MarkovMap(d2,[f1,f2],[0,0.5,1])
GC.acim(M2a)
@time rho2a = GC.acim(M1a)

pts = [points(d1);points(d2)]
@testapproxeq maxabs(rho1a(pts) - rho2a(pts)) 0.
@testapproxeq maxabs(rho1(pts) - rho2(pts)) 0.
