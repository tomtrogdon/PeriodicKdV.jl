# PeriodicKdV
Solve the IVP for periodic KdV and compute hyperelliptic solutions via dressing.

This respository contains the code for: [arxiv.org/abs/2205.00153](https://arxiv.org/abs/2205.00153)

An example usage follows.
```
gap_width = n -> isodd(n) ? exp(-n) : 0.1*n^5*exp(-n)
gap_start = n -> 4.0(n-1) + .4
g = 12
gaps = map(gap_start,1:g)
gaps = hcat(gaps, gaps + map(gap_width,1:g))
zs = gaps |> copy
zs[:,2] = 0*zs[:,2] .+ 1.0;
α1 = 0.1;
```
This code defines a sequence of intervals (`gaps`) and zeros `zs` that lie in these gaps along with their signs.  The variable `α1` gives the left most point the the spectrum of the associated hyperelliptic Schrodinger operator.  Then the code to do the precomputation to setup the associated KdV solution is given by:
```
S = HyperellipticSurface(gaps,zs,α1,300);
BA = BakerAkhiezerFunction(S,200.; tols = [1e-17,1e-10]);
u(x,t) = PeriodicKdV.KdV(BA,x,t,1e-8)
x = 0:.01:5
plot(x, map( x -> u(x,1.0), x) |> real)
```




[![CI](https://github.com/tomtrogdon/PeriodicKdV.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/tomtrogdon/PeriodicKdV.jl/actions/workflows/CI.yml) [![codecov](https://codecov.io/gh/tomtrogdon/PeriodicKdV.jl/branch/main/graph/badge.svg?token=JCU86U5O3J)](https://codecov.io/gh/tomtrogdon/PeriodicKdV.jl)

