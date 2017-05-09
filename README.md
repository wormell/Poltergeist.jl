# Poltergeist.jl 🏏👻

[![Build Status](https://travis-ci.org/wormell/Poltergeist.jl.svg?branch=master)](https://travis-ci.org/wormell/Poltergeist.jl)
[![Coverage Status](https://coveralls.io/repos/github/wormell/Poltergeist.jl/badge.svg?branch=master)](https://coveralls.io/github/wormell/Poltergeist.jl?branch=master)

Poltergeist is a package for quick, accurate and abstract approximation of statistical properties of one-dimensional chaotic dynamical systems. It is geared towards spectral methods (i.e. transfer operator-based approaches to dynamical systems) and uses spectral methods (e.g. Fourier and Chebyshev) for its numerical implementation. Poltergeist relies on and closely interfaces with the adaptive function approximation package [ApproxFun](https://github.com/ApproxFun/ApproxFun.jl).  

As an example your favourite Markov interval map and give it digital form:

```julia
using Poltergeist, ApproxFun
d = Segment(0,1)
fv = [x->2x+sin(2pi*x)/8,x->2-2x]
M = MarkovMap(fv,d)
M(0.25), M'(0.25)
```
<!---want to plot Markov Map--->

Similarly, take a circle map:
```julia
f_lift(x) = 2x + sin(2pi*x)/8pi
C = CircleMap(f_lift,PeriodicInterval(0,1))
```

Poltergeist's primary tool is a very efficient and usable computer representation of transfer operators. Calling ```Transfer``` on a ```MarkovMap``` automatically creates a calculable computer representation of the transfer operator in a Chebyshev basis. Transfer operators are instantiated as Operator types, and can do (numerically) all the kinds of things one expects from linear operators on function spaces:

```julia
L = Transfer(M)
f0 = Fun(x->sin(3pi*x),d) #ApproxFun function
f1 = L*f0
g = ((2I-L)\f0)'
``` 

Poltergeist has built-in commands for calculating many standard dynamical objects:

```julia
acim(C) # alias for acim(Transfer(C))
K = SolutionInv(L) # see Wormell, 2007
ρ = acim(K) # K caches QR factorisation of L
sigmasq = birkhoffvar(K,Fun(x->x^2,d))

using Plots
plot(ρ)
```
<!--- plot!(linearresponse(L,Fun(x->x*(1-x),d))) --->
<img src=https://github.com/johnwormell/Poltergeist.jl/raw/master/images/acim.png width=500 height=400>

One can even estimate eigenvalues:

```julia
eigvals(L,100)
```

Because the mathematical objects are given highly accurate Fourier and Chebyshev approximations, these commands give spectral accuracy: one can calculate acims to 15 digits of accuracy in a fraction of a second.


## Publications

This package is based on academic work. If you find this package useful in your work, please kindly cite as appropriate:

* J. P. Wormell (2017), Spectral collocation methods for transfer operators in uniformly expanding dynamics (preprint)

<!---
* S. Olver & A. Townsend (2014), A practical framework for infinite-dimensional linear algebra, Proceedings of the 1st First Workshop for High Performance Technical Computing in Dynamic Languages, 57–62
* J. P. Wormell (2017 in preparation), Fast numerical methods for intermittent systems 
--->

