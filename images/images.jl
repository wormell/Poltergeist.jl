using Poltergeist, Plots, ApproxFun; pyplot()
d = Interval(0,1);
M1(x) = 2x+sin(2pi*x)/6; M2(x) = 2-2x
M = MarkovMap([M1,M2],[0..0.5,0.5..1]);
L = Transfer(M)

scatter(eigvals(L,80),legend=true,grid=false,legend=false,
title = "Eigenvalues of transfer operator")
xlims!(-1.1,1.1)
ylims!(-0.5,0.5)
savefig("eigvals.pdf")
run(`sips -s format png eigvals.pdf --out eigvals.png`)

K = SolutionInv(L);
ρ = acim(K);
dρ = linearresponse(K,sinpi)
ϵ = 0.05
Mϵ = perturb(M,sinpi,ϵ)
plot(ρ,legend=true,grid=false,label="\$\\rho\_\{0\}\$",
title = "Linear response");
plot!(ρ+ϵ*dρ,label="\$\\rho\_\{\\epsilon\}\$ estimate")
plot!(acim(Mϵ), label="\$\\rho\_\{\\epsilon\}\$ actual")
xlabel!("\$x\$"); ylabel!("\$\\rho(x)\$")
# xlabel!("x"); ylabel!("ρ(x)")
ylims!(0.,1.8)
savefig("acim.pdf")
run(`sips -s format png acim.pdf --out acim.png`)
