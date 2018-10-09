# Pkg.installed()["ApproxFun"] != v"0.4.0+" && Pkg.checkout("ApproxFun","4bcc8f585361184342bb21780cc6be9893d99ce6")
using Poltergeist
using Base.Test
using ApproxFun, LinearAlgebra

f1(x)=2x+sin(2pi*x)/8pi; f2(x)=2x+sin(2pi*x)/8pi-1
f1d(x)=2+cos(2pi*x)/4; f2d = f1d

fv1(x) = x/2+sin(2pi*x)/8pi; fv2(x) = x/2+1/2+sin(2pi*x)/8pi
fv1d(x) = 1/2+cos(2pi*x)/4; fv2d = fv1d

#Periodic domain
println("Fourier tests 🌚🌞")
d1 = PeriodicInterval(0,1.)
M1b = CircleMap(fv1,d1,dir="rev",diff=fv1d)
L1b = Transfer(M1b)
@time L1b = Transfer(M1b)
K1b = SolutionInv(L1b)
acim(L1b)
@time ρ1b = acim(L1b)
acim(K1b)
@time ρ1b = acim(K1b)

acim(M1b)
@time ρ1b = acim(M1b)

M1f = CircleMap(f1,d1,diff=f1d)
acim(M1f)
@time ρ1f = acim(M1f)
println("Should all be ≤0.3s")

# Non-periodic domain
println("Chebyshev tests 🌳")
d2 = Segment(0..1.)
@test Poltergeist.coveringsegment([0..0.5,0.5..1]) == d2
M2b = MarkovMap([fv1,fv2],[0..0.5,0.5..1],dir=Reverse,diff=[fv1d,fv2d]);
acim(M2b)
@time ρ2b = acim(M2b)

M2ba = MarkovMap([fv1,fv2],[0..0.5,0.5..1],dir=Reverse); #autodiff comparison
acim(M2ba)
@time ρ2ba = acim(M2ba)

M2f = MarkovMap([f1,f2],[0..0.5,0.5..1],d2)
acim(M2f)
@time ρ2f = acim(M2f)
println("Should be ≤0.12s")

pts = [points(space(ρ1b),100);points(space(ρ2b),100)]
@test maximum(abs.(ρ1f.(pts) - ρ2f.(pts))) < 2000eps(1.)
@test maximum(abs.(ρ1b.(pts) - ρ2b.(pts))) < 2000eps(1.)
@test maximum(abs.(ρ2b.(pts) - ρ2ba.(pts))) < 2000eps(1.)

# # Transfer
# @test transfer(M1f,x->Fun(Fourier(d1),[0.,1.])(x),0.28531) == Poltergeist.transferfunction(0.28531,M1f,Poltergeist.BasisFun(Fourier(d1),2),Float64)
# @test transfer(M2f,exp,0.28531) ≈ (Transfer(M2f)*Fun(exp,Space(d2)))(0.28531)

println("Lanford map test")
lan = lanford()
K = SolutionInv(lan);
rho = acim(K);
l_exp = sum(Fun(x->log(abs(lan'(x))),0..1) * rho)
sigmasq_A = birkhoffvar(K,Fun(x->x^2,0..1))
L_lan = Transfer(lan)
K = SolutionInv(L_lan);
@time rho = acim(K);
@time l_exp = lyapunov(K)
@time l_exp2 = sum(Fun(x->log(abs(lan'(x))),0..1) * rho)
@time sigmasq_A = birkhoffvar(K,Fun(x->x^2,0..1))

@test l_exp ≈ 0.657661780006597677541582
@test l_exp2 ≈ 0.657661780006597677541582
@test sigmasq_A ≈ 0.360109486199160672898824

# modulomap test
println("Modulomap and examples test")
lan_lift(x) = 5x/2 - x^2/2
lan = modulomap(lan_lift,0..1);
@test Transfer(lan)[1:100,1:100] ≈ L_lan[1:100,1:100]
@test diag(Transfer(doubling(PeriodicInterval(6.,7.)))[1:10,1:10]) ≈ [1.;zeros(9)]
@test diag(Transfer(tupling(-4,0..4.))[1:10,1:10]) ≈ (-1/4).^(0:9)

# @test diag(Transfer(modulomap(x->1-x/5,0..1,dir=Reverse))[1:10,1:10]) .≈ (-0.2).^(0:9)

# Composing test
println("Composition test 🎼")
shiftmap = modulomap(x->5x+30,0..1,30..35)
lanshift = modulomap(x->5(x/5-6)/2-(x/5-6)^2/2,30..35,0..1)
@time doublelan = lan ∘ lanshift ∘ shiftmap
println("Should be ≤0.07s")
doubleK = SolutionInv(doublelan)
doublerho = acim(doubleK)
@time doublerho = acim(doublelan)
println("Should be ≤0.4s")
@test doublerho ≈ rho
@test lyapunov(doubleK) ≈ 2l_exp

lanpet = perturb(lan,sinpi,-0.1)∘inv(perturb(0..1,sinpi,-0.1))
@test lyapunov(lanpet) ≈ l_exp
@time lyapunov(lanpet)
println("Should be ≤0.01s")
# @time lanpet = acim(lanpet)
# @time lyapunov(lanpet)

# Eigvals test
println("Eigenvalue test")
c = 1/π
intervalmap = MarkovMap([x->sin(c*asin(x)),x->sin(c+(1-c)*asin(x))],[0..sin(c),sin(c)..sin(1.)],dir=Reverse)
eigs(intervalmap,100)
@time evs = eigvals(intervalmap,100)
println(length(evs))
@assert all(abs.(sort(evs,by=abs,rev=true)[1:5] - (c.^(1:5) + (1-c).^(1:5))).<1e-7)

# Correlation sums
println("Correlation sum test")
A1 = Fun(x->sin(sin(2pi*x)),d1)
A2 = Fun(x->sin(sin(2pi*x)),d2)
cs1f = correlationsum(M1f,A1)
@test maximum(abs.(cs1f.(pts)-correlationsum(M2f,A2).(pts))) .< 2000eps(1.)

println("Correlation function test")
A = Fun(x->x^2,0..1); B = Fun(sin,0..1)
cA,cB = covariancefunction(lan,A,B)
@test sum(cA)+sum(cB[2:end]) ≈ birkhoffcov(lan,A,B)
@time covariancefunction(lan,A,B)
  println("Should be <0.25s")
lancov = covariancefunction(lan,A)
@test lancov[1] + 2sum(lancov[2:end]) ≈ birkhoffvar(lan,A)
@time covariancefunction(lan,A)
covariancefunction(lan,A,100)
@time covariancefunction(lan,A,100)
  println("Should be <0.15s")

# Calling
println("Newton's method test ☏")
test_f = range(d2.a,stop=d2.b,length=20)[1:end-1] # map boundaries are dodgy because multivalued
test_x = [Poltergeist.mapinv(M2b,1,tf) for tf in test_f]
 @test M2b.(test_x) ≈ test_f
 @test M1b.(test_x) ≈ test_f
 @test M2b'.(test_x) ≈ M1b'.(test_x)

# #Inducing
# println("Inducing tests 🐴")
# M2bd = MarkovMap([fv1,fv2],[0..0.5,0.5..1],d2,dir=Reverse,diff=[fv1d,fv2d]);
# M2bi = induce(M2bd,1)
# # acim(M2bi)
# @time ρ2bi = acim(M2bi); println("Should be ≤4s")
# pts = points(space(ρ2bi),100)
#  normi = diff(cumsum(ρ2b).(∂(domain(M2bi))))[1]
# @test all(ρ2bi.(pts) .≈ ρ2b.(pts)/normi)# < 1000eps(1.)

# Time series
println("Time series tests")
NI = 10^6; NB = 10^3
@time ts = timeseries(M1f,NI,ρ1f)
println("Should be ≤2s")
# println(ts[1:11])
# println(ts[end-10:end])
@test abs(sum(sin.(sinpi.(2ts)))/NI #- sum(ρ1f*A1)  # which is 0
      )< (4sum(cs1f*A1)+200eps(1.))/sqrt(NI)

@time cts = timehist(M2f,NI,NB,ρ2f)
@test abs(sum(sin.(sin.(2pi*cts[1][1:end-1])).*cts[2])/NI - sum(ρ2f*A2))< 1/NB+(4sum(cs1f*A1)+200eps(1.))/sqrt(NI)
println("Should be ≤2s")

#TODO: fix
# # Intermittent maps
# println("Intermittent tests")
# for α in [0.22,1.3523]
#   println("α = $α")
#   @time b = NeutralBranch(x->1+2^α*x,x->2^α,α,0.6/2^α,Interval(0,0.5),Interval(0,1))
#   # @time b = NeutralBranch(x->1+2^α*x,x->2^α,α,0.6/2^α,Interval(0,0.5),Interval(0,1))
#   b2 = branch(x->(x+1)/2,x->0.5,Segment(0.5,1.),Segment(0.,1.),dir="rev")
#   Mint = MarkovMap(Segment(0.,1),Segment(0.,1),[b,b2])
#
#   Mint_I = induce(Mint,1)
#   ρint = acim(Transfer(Mint_I))
#   L = Transfer(Mint_I)
#   @time ρint = acim(L)
#   println("Size of transfer operator: $(L.datasize)")
#   # TODO: FullAcim
#   # pts = points(domainspace(L),40)
#   # @time pfull =
# end

# TODO: put back in on new ApproxFun release (>0.9.0)
# 2D tests - in testing
# println("2D tests")
# using StaticArrays
# standardmap_inv_lift(x::SVector) = SVector(x[1] - 0.1*sin(x[2] - x[1]),x[2]-x[1]);
# standardmap_inv_diff(x::SVector) = SMatrix{2,2}(1,0,0,1); # As only determinant is important...
# dom = PeriodicInterval()^2
#  # binv = branch(standardmap_inv_lift,standardmap_inv_diff,dom,dom,dir="rev"); # deprecated
# binv= branch(standardmap_inv_lift,dom,dom,standardmap_inv_diff,dir=Reverse)
# standardmap = MarkovMap([binv],dom,dom)
# L_standard = Transfer(standardmap)
# ApproxFun.resizedata!(L_standard,:,2)

println("")
println("😎")
