NonOpt
======

This is a RELEASE CANDIDATE of NonOpt (version 1.0rc3), meaning that it may still contain significant bugs that have yet to emerge through testing.  If you encounter any significant bug(s), then please contact [Frank E. Curtis](mailto:frank.e.curtis@gmail.com).

Overview
--------

NonOpt (Nonlinear/nonconvex/nonsmooth Optimizer) is a software package for solving minimization problems.  It is designed to locate a stationary point (ideally, a minimizer) of

```
min     f(x)
x ∈ Rⁿ
```

where ```f : Rⁿ -> R``` is locally Lipschitz and continuously differentiable over a full-measure subset of ```Rⁿ```.  The function ```f``` can be nonlinear, nonconvex, and/or nonsmooth.

NonOpt is written in C++ and is released under the MIT License.  The main author is [Frank E. Curtis](http://coral.ise.lehigh.edu/frankecurtis/).  For a list of all contributors, please see the [AUTHORS file](NonOpt/AUTHORS).

Compiling NonOpt requires BLAS and LAPACK.  The code for these packages is not provided in this repository, which are available under different conditions and licenses than those for NonOpt.

Please visit the [NonOpt homepage](http://coral.ise.lehigh.edu/frankecurtis/nonopt/).

Citing NonOpt
-------------

NonOpt is provided free of charge so that it might be useful to others.  Please send e-mail to [Frank E. Curtis](mailto:frank.e.curtis@gmail.com) with success stories or other feedback.  If you use NonOpt in your research, then, for now, please cite the following papers (or at least the first one):

- F. E. Curtis and M. Li. "Gradient Sampling Methods with Inexact Subproblem Solutions and Gradient Aggregation." [arXiv:2005.07822](https://arxiv.org/abs/2005.07822), 2020.
- F. E. Curtis, D. P. Robinson, and B. Zhou. "A Self-Correcting Variable-Metric Algorithm Framework for Nonsmooth Optimization." IMA Journal of Numerical Analysis, 10.1093/imanum/drz008, 2019.
- F. E. Curtis and X. Que. "A Quasi-Newton Algorithm for Nonconvex, Nonsmooth Optimization with Global Convergence Guarantees." Mathematical Programming Computation, 7(4):399–428, 2015.
- F. E. Curtis and X. Que. "An Adaptive Gradient Sampling Algorithm for Nonsmooth Optimization." Optimization Methods and Software, 28(6):1302–1324, 2013.
