using Poltergeist, Plots, ApproxFun; pyplot()
d = Segment(0,1);
M1(x) = 2x+sin(2pi*x)/6; M2(x) = 2-2x
M = MarkovMap([M1,M2],[0..0.5,0.5..1]);
K = SolutionInv(M);
ρ = acim(K);
dρ = linearresponse(K,sinpi)
ϵ = 0.05
ϵpert(x) = x + ϵ*sinpi(x)
Mϵ = MarkovMap([ϵpert∘M1,ϵpert∘M2],[0..0.5,0.5..1])
plot(ρ,legend=true,grid=false,label="\$\\rho\_\{0\}\$",
title = "Linear response");
plot!(ρ+ϵ*dρ,label="\$\\rho\_\{\\epsilon\}\$ estimate")
plot!(acim(Mϵ), label="\$\\rho\_\{\\epsilon\}\$ actual")
xlabel!("\$x\$"); ylabel!("\$\\rho(x)\$")
# xlabel!("x"); ylabel!("ρ(x)")
ylims!(0.,1.8)
savefig("acim.pdf")
println("First image done")
