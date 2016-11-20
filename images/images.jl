using GhostCoop, Plots, ApproxFun; pyplot()
d = Interval(0,1)
M = MarkovMap(d,[x->2x+sin(2pi*x)/6,x->2-2x],[0,1/2,1])
L = Transfer(M)
ρ = acim(L) # acim(M) also works but L caches itself
plot(ρ,legend=false,grid=false);
xlabel!("\$x\$"); ylabel!("\$\\rho(x)\$")
ylims!(0.,1.5)
PyPlot.savefig("acim.png",dpi=300)
println("First image done")