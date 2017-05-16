# Poltergeist.jl 

[![Build Status](https://travis-ci.org/wormell/Poltergeist.jl.svg?branch=master)](https://travis-ci.org/wormell/Poltergeist.jl)
[![Coverage Status](https://coveralls.io/repos/github/wormell/Poltergeist.jl/badge.svg?branch=master)](https://coveralls.io/github/wormell/Poltergeist.jl?branch=master)

Poltergeist is a package for quick, accurate and abstract approximation of statistical properties of one-dimensional chaotic dynamical systems. 

It treats chaotic systems through the framework of spectral methods (i.e. transfer operator-based approaches to dynamical systems) and is numerically implemented via spectral methods (e.g. Fourier and Chebyshev). For the latter, Poltergeist relies on and closely interfaces with the adaptive function approximation package [ApproxFun](https://github.com/ApproxFun/ApproxFun.jl).  

The use of highly accurate Fourier and Chebyshev approximations, means that we have spectral convergence rates: one can calculate acims to 15 digits of accuracy in a fraction of a second.

As an example, take your favourite Markov interval map and give it digital form:

```julia
using Poltergeist, ApproxFun
d = Segment(0,1)
f1 = 2x+sin(2pi*x)/8; f2(x) = 2-2x
f = MarkovMap([f1,f2],[0..0.5,0.5..1])
f(0.25), f'(0.25)
```
<!---want to plot Markov Map--->

Similarly, take a circle map, or maps defined by modulo or inverse:

```julia
c = CircleMap(x->4x + sin(2pi*x)/2pi,PeriodicInterval(0,1))
lanford = modulomap(x->2x+x*(1-x)/2,0..1)
doubling = MarkovMap([x->x/2,x->(x+1)/2],0..1,dir=Reverse)
```

<!---For better performance, use generic (vs anonymous) functions and (if using a complicated function) supply a derivative:
```julia
complicatedfun(x) = 3x+sum(2^(-33/8)m)
--->

Calling ```Transfer``` on an ```AbstractMarkovMap``` type automatically creates an ApproxFun ```Operator```, with which you can do (numerically) all the kinds of things one expects from linear operators on function spaces:

```julia
L = Transfer(M)
f0 = Fun(x->sin(3pi*x),d) #ApproxFun function
f1 = L*f0
g = ((2I-L)\f0)'
eigvals(L,30)
det(I-4L) # Fredholm determinant
``` 

In particular, you can solve for many statistical properties, many of which Poltergeist has built-in commands for. Most of these commands allow you to use the ```MarkovMap``` directly (bad, zero caching between uses), transfer operator (caches transfer operator entries, usually the slowest step), or the ```SolutionInv``` operator (caches QR factorisation as well).

```julia
K = SolutionInv(L)
ρ = acim(K)
@test ρ == K\Fun(one,d)
birkhoffvar(K,Fun(x->x^2,d))
birkhoffcov(K,Fun(x->x^2,d),Fun(identity,d))
dρ = linearresponse(K,Fun(sinpi,d))

using Plots
plot(ρ)
ε = 0.05
plot!(ρ + ε*dρ)
pertε(x) = x + ε*sinpi(x)
plot!(acim(MarkovMap([pertε∘f1,pertε∘f2],[0..0.5,0.5..1])))
```
<!--- TODO: plot!(linearresponse(L,Fun(x->x*(1-x),d))) --->
<img src=https://github.com/johnwormell/Poltergeist.jl/raw/master/images/acim.pdf width=500 height=400>

## Publications

This package is based on academic work. If you find this package useful in your work, please kindly cite as appropriate:

* J. P. Wormell (2017), Spectral Galerkin methods for transfer operators in uniformly expanding dynamics ([preprint](https://arxiv.org/abs/1705.04431))
<!---* *features of ApproxFun*: S. Olver & A. Townsend (2014), A practical framework for infinite-dimensional linear algebra, Proceedings of the 1st First Workshop for High Performance Technical Computing in Dynamic Languages, 57–62 ([article](https://doi.org/10.1109/HPTCDL.2014.10), [preprint](https://arxiv.org/abs/1409.5529)) --->
* Literature [on ApproxFun](https://github.com/JuliaApproximation/ApproxFun.jl#references).


<!---
* J. P. Wormell (2017 in preparation), Fast numerical methods for intermittent systems 


_____________
* A map f on [-1,1] is C-expanding if acos◦f◦cos is uniformly expanding. Most uniformly expanding maps are C-expanding, or if not there is always a conjugate or iterate that is.

--->