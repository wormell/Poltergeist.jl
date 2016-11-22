Pkg.installed()["ApproxFun"] != v"0.4.0+" && Pkg.checkout("ApproxFun","development")
using Poltergeist
using Base.Test
using ApproxFun

f1(x)=2x+sin(2pi*x)/8pi; f2(x)=2x+sin(2pi*x)/8pi-1
f1d(x)=2+cos(2pi*x)/4; f2d = f1d

fv1(x) = x/2+sin(2pi*x)/8pi; fv2(x) = x/2+1/2+sin(2pi*x)/8pi
fv1d(x) = 1/2+cos(2pi*x)/4; fv2d = fv1d

#Periodic domain
println("Fourier tests")
d1 = PeriodicInterval([0,1])
M1 = (MarkovMap(d1,[fv1,fv2],[fv1d,fv2d],[0,0.5,1],"rev"));
acim(M1)
@time rho1b = acim(M1)

M1a = MarkovMap(d1,[f1,f2],[f1d,f2d],[0,0.5,1])
acim(M1a)
@time rho1f = acim(M1a)

# Non-periodic domain
println("Chebyshev tests")
d2 = Interval([0,1])
M2b = MarkovMap(d2,[fv1,fv2],[fv1d,fv2d],[0,0.5,1],"rev");
acim(M2b)
@time rho2b = acim(M2b)

M2ba = MarkovMap(d2,[fv1,fv2],[0,0.5,1],"rev"); #autodiff comparison
acim(M2ba)
Profile.clear()
@time rho2ba = acim(M2ba)



M2f = MarkovMap(d2,[f1,f2],[0,0.5,1])
acim(M2f)
@time rho2f = acim(M2f)

pts = [points(space(rho1b),100);points(space(rho2b),100)]
@test maxabs(rho1f(pts) - rho2f(pts)) < 200eps(1.)
@test maxabs(rho1b(pts) - rho2b(pts)) < 200eps(1.)
@test maxabs(rho2b(pts) - rho2ba(pts)) < 200eps(1.)


#Inducing
println("Inducing tests")
M2bd = MarkovMap(d2,[fv1,fv2],[fv1d,fv2d],[0,0.5,1],"rev");
M2bi = induce(M2bd,1)
# acim(M2bi)
@time rho2bi = acim(M2bi)
pts = points(space(rho2bi),100)
normi = diff(cumsum(rho2b)(∂(domain(M2bi))))[1]
@test maxabs(rho2bi(pts) - rho2b(pts)/normi) < 200eps(1.)

