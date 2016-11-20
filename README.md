# GhostCoop 👻✊

[![Build Status](https://travis-ci.org/johnwormell/GhostCoop.jl.svg?branch=master)](https://travis-ci.org/johnwormell/GhostCoop.jl)
[![Coverage Status](https://coveralls.io/repos/github/johnwormell/GhostCoop.jl/badge.svg?branch=master)](https://coveralls.io/github/johnwormell/GhostCoop.jl?branch=master)

GhostCoop is a package for quick and accurate approximation of one-dimensional chaotic dynamical systems. It is geared towards spectral methods* and uses spectral methods^† for its numerical implementation.


 <!---computes transfer operators of one-dimensional chaotic systems in spectral bases. This enables one to find statistical properties of dynamical systems quickly, reliably and abstractly.
--->

Take your favourite Markov uniformly expanding dynamical system and create a numerical representation of it:

```julia
using GhostCoop, ApproxFun
d = Interval(0,1)
M = MarkovMap(d,[x->2x+sin(2pi*x)/8,x->2-2x],[0,1/2,1])
M(0.25), M'(0.25)
```
<!---want to plot Markov Map--->

GhostCoop closely interfaces with the adaptive function approximation package [ApproxFun](https://github.com/ApproxFun/ApproxFun.jl). Calling ```Transfer``` on a ```MarkovMap``` automatically creates a calculable computer representation of the transfer operator in a Chebyshev basis. Transfer operators are instantiated as Operator types, and can do (numerically) all the kinds of things one expects from linear operators on function spaces:

```julia
L = Transfer(M)
f0 = Fun(x->sin(3pi*x),d) #ApproxFun function
f1 = L*f0
g = ((2I-L)\f0)'
sum(f1*Fun(sin,d)) # sum is total integral
sum(f1*Fun(sin,d)) - sum(f0*Fun(x->sin(M(x)),d)) #definition of transfer operator as adjoint
``` 

GhostCoop has built-in commands for calculating many standard dynamical objects:

```julia
ρ = acim(L) # acim(M) also works but L is already cached
norm(L*ρ-ρ)
using Plots
plot(ρ)
```
<!--- plot!(linearresponse(L,Fun(x->x*(1-x),d))) --->
<img src=https://github.com/johnwormell/GhostCoop.jl/raw/master/images/acim.png width=500 height=400>

These have spectral accuracy: one can calculate acims to 15 digits of accuracy in a fraction of a second:

```julia
doubling = MarkovMap(Interval(0,1),[x->2x,x->2x-1],[0,0.5,1])
norm(acim(doubling)-1) #doubling map
```


## Publications

This package is based on academic work that is being prepared for submission. If you find this package useful in your work, please kindly cite the following papers as appropriate. Please also send me an email (<j.wormell@maths.usyd.edu.au>), I'd love to hear about it! 

* J. P. Wormell (2016 in preparation), Spectral collocation methods for transfer operators in uniformly expanding dynamics
* S. Olver & A. Townsend (2014), A practical framework for infinite-dimensional linear algebra, Proceedings of the 1st First Workshop for High Performance Technical Computing in Dynamic Languages, 57–62

<!--- J. P. Wormell (2017 in preparation), Fast numerical methods for intermittent systems --->



________________
<sub>*  Transfer operator-based approaches to dynamical systems</sub>

<sub>† Fourier and Chebyshev</sub>
