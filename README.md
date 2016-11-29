# Poltergeist.jl 🏏👻

[![Build Status](https://travis-ci.org/wormell/Poltergeist.jl.svg?branch=master)](https://travis-ci.org/wormell/Poltergeist.jl)
[![Coverage Status](https://coveralls.io/repos/github/wormell/Poltergeist.jl/badge.svg?branch=master)](https://coveralls.io/github/wormell/Poltergeist.jl?branch=master)

Poltergeist is a package for quick and accurate approximation of one-dimensional chaotic dynamical systems. It is geared towards spectral methods* and uses spectral methods^† for its numerical implementation. Poltergeist relies on and closely interfaces with the adaptive function approximation package [ApproxFun](https://github.com/ApproxFun/ApproxFun.jl).  


 <!---computes transfer operators of one-dimensional chaotic systems in spectral bases. This enables one to find statistical properties of dynamical systems quickly, reliably and abstractly.
--->

<!---Take your favourite Markov uniformly expanding dynamical system and give it digital form:

```julia
using Poltergeist, ApproxFun
d = Interval(0,1)
fv = [x->2x+sin(2pi*x)/8,x->2-2x]
M = MarkovMap(d,fv,[0,1/2,1])
M(0.25), M'(0.25)
```
<>---want to plot Markov Map---/>

Primarily, Poltergeist provides a very efficient and usable computer representation of transfer operators. Calling ```Transfer``` on a ```MarkovMap``` automatically creates a calculable computer representation of the transfer operator in a Chebyshev basis. Transfer operators are instantiated as Operator types, and can do (numerically) all the kinds of things one expects from linear operators on function spaces:

```julia
L = Transfer(M)
f0 = Fun(x->sin(3pi*x),d) #ApproxFun function
f1 = L*f0
g = ((2I-L)\f0)'
sum(f1*Fun(sin,d)) # sum is total integral
sum(f1*Fun(sin,d)) - sum(f0*Fun(x->sin(M(x)),d)) #definition of transfer operator as adjoint
``` 

Poltergeist has built-in commands for calculating many standard dynamical objects:

```julia
ρ = acim(L) # acim(M) also works but L is already cached
norm(L*ρ-ρ) # 2.3043432e-16
using Plots
plot(ρ)
```
<!--- plot!(linearresponse(L,Fun(x->x*(1-x),d))) ---/>
<img src=https://github.com/johnwormell/GhostCoop.jl/raw/master/images/acim.png width=500 height=400>

Because the mathematical objects are given highly accurate Fourier and Chebyshev approximations, these commands give spectral accuracy: one can calculate acims to 15 digits of accuracy in a fraction of a second.

```julia
fdv = [x->2+2pi*cos(2pi*x)/8,x->-2.]
Md = MarkovMap(d,fv,fdv,[0,1/2,1]) # providing a derivative speeds things up a lot
@time acim(M) # 0.686850 seconds
```
<!--- ```julia
doubling = MarkovMap(Interval(0,1),[x->2x,x->2x-1],[0,0.5,1])  #doubling map
norm(acim(doubling)-1)
```
--->

<!--- induced maps ---/>

## Publications

This package is based on academic work that is being prepared for submission. If you find this package useful in your work, please kindly cite the following papers as appropriate. Seeing as some papers are in preparation, please also send me an email to check on progress (<j.wormell@maths.usyd.edu.au>)!

* J. P. Wormell (2016 in preparation), Spectral collocation methods for transfer operators in uniformly expanding dynamics
* S. Olver & A. Townsend (2014), A practical framework for infinite-dimensional linear algebra, Proceedings of the 1st First Workshop for High Performance Technical Computing in Dynamic Languages, 57–62

<!--- J. P. Wormell (2017 in preparation), Fast numerical methods for intermittent systems --->



________________
<sub>*  Transfer operator-based approaches to dynamical systems</sub>

<sub>† Fourier and Chebyshev</sub>
