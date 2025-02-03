NonOpt
======

This is a version 2.0 of NonOpt.  If you encounter any significant bug(s), then please send e-mail to [Frank E. Curtis](mailto:frank.e.curtis@gmail.com).

Overview
--------

NonOpt (Nonlinear/nonconvex/nonsmooth Optimizer) is a software package for solving minimization problems.  It is designed to locate a stationary point (ideally, a minimizer) of

```
min     f(x)
x ∈ Rⁿ
```

where ```f : Rⁿ -> R``` is locally Lipschitz over ```Rⁿ``` and continuously differentiable over a full-measure subset of ```Rⁿ```.  The function ```f``` can be nonlinear, nonconvex, and/or nonsmooth.

NonOpt is written in C++ and is released under the MIT License.  The main author is [Frank E. Curtis](http://coral.ise.lehigh.edu/frankecurtis/).  For a list of all contributors, please see the [AUTHORS file](NonOpt/AUTHORS).

Compiling NonOpt requires BLAS and LAPACK.  The code for these packages is not provided in this repository and must be obtained separately by a user.

Please visit the [NonOpt homepage](http://coral.ise.lehigh.edu/frankecurtis/nonopt/).

Citing NonOpt
-------------

NonOpt is provided free of charge so that it might be useful to others.  Please send e-mail to [Frank E. Curtis](mailto:frank.e.curtis@gmail.com) with success stories or other feedback.  If you use NonOpt in your research, then please cite the following papers:

- Frank E. Curtis and Lara Zebiane. "NonOpt: Nonconvex, Nonsmooth Optimizer." 2025.
- Frank E. Curtis and Minhan Li. "Gradient Sampling Methods with Inexact Subproblem Solutions and Gradient Aggregation." INFORMS Journal on Optimization, 4(4):426-445, 2022.
- Frank E. Curtis, Daniel P. Robinson, and Baoyu Zhou. "A Self-Correcting Variable-Metric Algorithm Framework for Nonsmooth Optimization." IMA Journal of Numerical Analysis, 40(2):1154-1187, 2020.
- Frank E. Curtis and Xiaocun Que. "A Quasi-Newton Algorithm for Nonconvex, Nonsmooth Optimization with Global Convergence Guarantees." Mathematical Programming Computation, 7(4):399–428, 2015.
- Frank E. Curtis and Xiaocun Que. "An Adaptive Gradient Sampling Algorithm for Nonsmooth Optimization." Optimization Methods and Software, 28(6):1302–1324, 2013.
