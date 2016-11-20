using GhostCoop
using Base.Test

Pkg.checkout("ApproxFun","development")
using ApproxFun

# write your own tests here

# Expanding tests

f1(x)=2x+sin(2pi*x)/8pi; f2(x)=2x+sin(2pi*x)/8pi-1
fv1(x) = x/2+sin(2pi*x)/8pi; fv2(x) = x/2+1/2+sin(2pi*x)/8pi

# Periodic domain
println("Fourier tests")
d1 = PeriodicInterval([0,1])
M1 = (MarkovMap(d1,[fv1,fv2],[0,0.5,1],"rev"));
acim(M1)
@time rho1 = acim(M1)

M1a = MarkovMap(d1,[f1,f2],[0,0.5,1])
acim(M1a)
@time rho1a = acim(M1a)

# Non-periodic domain
println("Chebyshev tests")
d2 = Interval([0,1])
M2 = MarkovMap(d2,[fv1,fv2],[0,0.5,1],"rev");
acim(M2)
@time rho2 = acim(M2)

M2a = MarkovMap(d2,[f1,f2],[0,0.5,1])
acim(M2a)
@time rho2a = acim(M2a)

pts = [points(space(rho1),100);points(space(rho2),100)]
@test maxabs(rho1a(pts) - rho2a(pts)) < 200eps(1.)
@test maxabs(rho1(pts) - rho2(pts)) < 200eps(1.)
