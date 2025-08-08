# IncSignFeas.jl

[![Build Status](https://github.com/BP-J/IncSignFeas.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BP-J/IncSignFeas.jl/actions/workflows/CI.yml?query=branch%3Amain)

`julia` code aimed at identifying the chambers of an hyperplane arrangement (knowing the sign vectors of the chambers), currently being in work to be added as a package.

This package uses (in most algorithms), a tree structure introduced by [Rada and Černý](https://epubs.siam.org/doi/10.1137/15M1027930). 
Several variants of this tree algorithms are proposed, some relying on linear optimization and some on stem vectors (a variant of signed circuits of matroids).
For additional details, see [here](https://inria.hal.science/hal-05002249). 

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
Then, `options` defining which algorithm will be used must be defined, for instance with `options = options_from_algo(3, false)`.
Finally, the algorithm is called with `info = isf([V ; t'], options)`. 