# IncSignFeas.jl

[![Build Status](https://github.com/BP-J/IncSignFeas.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BP-J/IncSignFeas.jl/actions/workflows/CI.yml?query=branch%3Amain)

IncSignFeas (abbreviated as `isf`) is a `julia` code aimed at identifying the chambers of an hyperplane arrangement (knowing the sign vectors of the chambers).

This package uses a tree structure introduced by [Rada and Černý](https://epubs.siam.org/doi/10.1137/15M1027930). 
Several variants of this tree algorithms are proposed, some relying on linear optimization and some on stem vectors (a variant of signed circuits of matroids).
For additional details, see [here](https://inria.hal.science/hal-05002249) and [here](https://bp-j.github.io/research/documentation_ISFjl.pdf). 
To output the tables relevant in these documents, see the very end of this readme. 

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

The main tools used require the `julia` packages [HiGHS](https://github.com/jump-dev/HiGHS.jl) ([Gurobi](https://github.com/jump-dev/Gurobi.jl) was used by the authors, though any linear optimization solver such as [GLPK](https://github.com/jump-dev/GLPK.jl) works too), [JuMP](https://github.com/jump-dev/JuMP.jl), though some algorithms can work without these. 

## Main principle

An arrangement is defined by a matrix `V` of size **n x p** and a vector `t` of size **p**, and must be grouped as `[V ; t']`.
The arrangement contains **p** hyperplanes defined by **{x : V[:,i]'x = t[i]}**.

The code IncSignFeas (often abbreviated ISF or isf) contains a few different parameters defined as `options`, 
the main one being the algorithm that will be used. We give a few examples. The instructions
```
Vt = [1 0 1 1 ; 0 1 1 -1 ; 0 0 1 0];
options = options_from_algo(13, false); # overall best algorithm
info = isf(Vt, options);
info.ns
info.s
```
generate an arrangement, the options (`false` indicates the arrangement is not central), and launches the code. 
Then, the number of sign vectors and the set of sign vectors are displayed, which should be 10 and 
`[[+1,+1,+1,+1], [+1,+1,+1,-1], [+1,+1,-1,+1], [+1,+1,-1,-1], [+1,-1,-1,+1], [+1,-1,+1,+1], [-1,+1,-1,-1], [-1,+1,+1,-1], [-1,-1,-1,+1], [-1,-1,-1,-1]]`.

The main algorithms are obtained with the instructions (then use `info = isf(Vt, options);` for each)
```
options = options_from_algo(5, false);
options = options_from_algo(0, false); # the initial tree algorithm of Rada and Cerny
options = options_from_algo(3, false);
options = options_from_algo(11, false);
options = options_from_algo(7, false);
options = options_from_algo(15, false);
```
which change some parameters of `info` and mostly the order in which `info.s` is displayed.

To treat linear arrangements, input `V` (without a line of zeros for `t`) and replace `false` to `true` when defining the options. 

## Output "result" tables

In the folder test, the isf_benchmark_launches_aff.jl file is the main file on which the experiments were run. 
The file isf_tables.jl contains a code that produces tables very similar to those found at the end of the associated papers. 
After using `include("isf_tables.jl")`, the following command should produce the corresponding tables. 
```
table_72, table_A1, table_A2, table_A3 = tables_affines(list_short) # only a third of the instances
table_72, table_A1, table_A2, table_A3 = tables_affines(list_total) # quite long, all the instances
```
Those commands are intentionally quickened in the sense an algorithm on an instance is only ran 3 times (not more which can be relevant for faster executions to reduce variance). 