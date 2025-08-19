# IncSignFeas.jl

[![Build Status](https://github.com/BP-J/IncSignFeas.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BP-J/IncSignFeas.jl/actions/workflows/CI.yml?query=branch%3Amain)

`julia` code aimed at identifying the chambers of an hyperplane arrangement (knowing the sign vectors of the chambers), currently being in work to be added as a package.

This package uses (in most algorithms), a tree structure introduced by [Rada and Černý](https://epubs.siam.org/doi/10.1137/15M1027930). 
Several variants of this tree algorithms are proposed, some relying on linear optimization and some on stem vectors (a variant of signed circuits of matroids).
For additional details, see [here](https://inria.hal.science/hal-05002249) and [here](https://bp-j.github.io/research/documentation_ISFjl_1908.pdf). 

## Package presentation

The package can be loaded in `julia` by

```
]
add IncSignfeas
using IncSignFeas
```
or 
```
using Pkg; Pkg.add("IncSignFeas")
using IncSignFeas
```

The main tools used require the `julia` packages [Gurobi](https://github.com/jump-dev/Gurobi.jl) (or any linear optimization solver such as [GLPK](https://github.com/jump-dev/GLPK.jl) or [HiGHS](https://github.com/jump-dev/HiGHS.jl) which are license-free and available in `julia`), [JuMP](https://github.com/jump-dev/JuMP.jl), though some algorithms can work without these. 

## Main principle

An arrangement is defined by a matrix **V** of size **n x p** and a vector **t** of size **p**, and must be grouped as `[V ; t']`.
The arrangement contains **p** hyperplanes defined by **{x : V[:,i]'x = t[i]}**.

The code IncSignFeas (often abbreviated ISF or isf) contains a few different parameters defined as `options`, 
the main one being the algorithm that will be used. A test example can be as follows:
```
Vt = [1 0 1 1 ; 0 1 1 -1 ; 0 0 1 0]
options = options_from_algo(0, false)
info = isf(Vt, options)
info.s
```
which generates an arrangement, the options (`false` indicates the arrangement is not central), and launches the code. 
Then, the set of sign vectors is displayed, which should be 
`[[+1,+1,+1,+1], [+1,+1,+1,-1], [+1,+1,-1,+1], [+1,+1,-1,-1], [+1,-1,-1,+1], [+1,-1,+1,+1], [-1,+1,-1,-1], [-1,+1,+1,-1], [-1,-1,-1,+1], [-1,-1,-1,-1]]`.
