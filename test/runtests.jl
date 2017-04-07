# Pkg.installed()["ApproxFun"] != v"0.4.0+" && Pkg.checkout("ApproxFun","4bcc8f585361184342bb21780cc6be9893d99ce6")

using Poltergeist
using Base.Test, Compat
using ApproxFun

f1(x)=2x+sin(2pi*x)/8pi; f2(x)=2x+sin(2pi*x)/8pi-1
f1d(x)=2+cos(2pi*x)/4; f2d = f1d

fv1(x) = x/2+sin(2pi*x)/8pi; fv2(x) = x/2+1/2+sin(2pi*x)/8pi
fv1d(x) = 1/2+cos(2pi*x)/4; fv2d = fv1d

#Periodic domain
## TODO: put back in with CircleMap
# println("Fourier tests")
d1 = PeriodicInterval(0,1.)
# M1b = (MarkovMap([fv1,fv2],[fv1d,fv2d],[0..0.5,0.5..1],d1,"rev"));
# acim(M1b)
# @time ρ1b = acim(M1b)
#
# M1f = MarkovMap([f1,f2],[f1d,f2d],[0..0.5,0.5..1],d1)
# acim(M1f)
# @time ρ1f = acim(M1f)
#
# Non-periodic domain
println("Chebyshev tests")
d2 = Segment(d1)
M2b = MarkovMap([fv1,fv2],[fv1d,fv2d],[0..0.5,0.5..1],d2,dir="rev");
acim(M2b)
@time ρ2b = acim(M2b)

M2ba = MarkovMap([fv1,fv2],[0..0.5,0.5..1],d2,dir="rev"); #autodiff comparison
acim(M2ba)
@time ρ2ba = acim(M2ba)



M2f = MarkovMap([f1,f2],[0..0.5,0.5..1],d2)
acim(M2f)
@time ρ2f = acim(M2f)

## TODO: put back in with CircleMap
# pts = [points(space(ρ1b),100);points(space(ρ2b),100)]
# @test maxabs(ρ1f(pts) - ρ2f(pts)) < 200eps(1.)
# @test maxabs(ρ1b(pts) - ρ2b(pts)) < 200eps(1.)
# @test maxabs(ρ2b(pts) - ρ2ba(pts)) < 200eps(1.)
#
# # # Transfer
# # @test transfer(M1f,x->Fun(Fourier(d1),[0.,1.])(x),0.28531) == Poltergeist.transferfunction(0.28531,M1f,Poltergeist.BasisFun(Fourier(d1),2),Float64)
# # @test_approx_eq transfer(M2f,exp,0.28531) (Transfer(M2f)*Fun(exp,Space(d2)))(0.28531)
#
#
# # Correlation sums
# A1 = Fun(x->sin(sin(2pi*x)),d1)
# A2 = Fun(x->sin(sin(2pi*x)),d2)
# cs1f = correlationsum(M1f,A1)
# @test maxabs(cs1f(pts)-correlationsum(M2f,A2)(pts)) < 20000eps(1.)

#Inducing
println("Inducing tests")
M2bd = MarkovMap([fv1,fv2],[fv1d,fv2d],[0..0.5,0.5..1],d2,dir="rev");
M2bi = induce(M2bd,1)
# acim(M2bi)
@time ρ2bi = acim(M2bi)
pts = points(space(ρ2bi),100)
@compat normi = diff(cumsum(ρ2b).(∂(domain(M2bi))))[1]
@test maxabs(ρ2bi(pts) - ρ2b(pts)/normi) < 200eps(1.)

## TODO: put back in with CircleMap
# # Time series
# println("Time series tests")
# NI = 10^6; NB = 10^3
# @time ts = timeseries(M1f,NI,ρ1f)
# @test abs(sum(sin(sin(2pi*ts)))/NI - sum(ρ1f*A1))< (4sum(cs1f*A1)+200eps(1.))/sqrt(NI)
#
# @time cts = timehist(M2f,NI,NB,ρ2f)
# @test abs(sum(sin(sin(2pi*cts[1][1:end-1])).*cts[2])/NI - sum(ρ2f*A2))< 1/NB+(4sum(cs1f*A1)+200eps(1.))/sqrt(NI)

# Intermittent maps
println("Intermittent tests")
for α in [0.22,1.3523]
  println("α = $α")
  @time b = NeutralBranch(x->1+2^α*x,x->2^α,α,0.6/2^α,Interval(0,0.5),Interval(0,1))
  # @time b = NeutralBranch(x->1+2^α*x,x->2^α,α,0.6/2^α,Interval(0,0.5),Interval(0,1))
  b2 = branch(x->(x+1)/2,x->0.5,Segment(0.5,1.),Segment(0.,1.),dir="rev")
  Mint = MarkovMap(Segment(0.,1),Segment(0.,1),[b,b2])

  Mint_I = induce(Mint,1)
  ρint = acim(Transfer(Mint_I))
  L = Transfer(Mint_I)
  @time ρint = acim(L)
  println("Size of transfer operator: $(L.datasize)")
  # TODO: FullAcim
  # pts = points(domainspace(L),40)
  # @time pfull =
end
