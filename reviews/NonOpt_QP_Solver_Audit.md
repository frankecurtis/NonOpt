# NonOpt — QP Solver Audit & Build/CI Deliverables

**Repository:** `frankecurtis/NonOpt` (v2.1) · **Date:** July 12, 2026
**Files audited:** `NonOptQPSolverDualActiveSet.{hpp,cpp}` (2,426 LOC), `NonOptQPSolverInteriorPoint.{hpp,cpp}` (1,335 LOC), `NonOptQPSolver.hpp`.
**Companion deliverables:** a CMake build (`nonopt-build/`) and a GitHub Actions workflow (`nonopt-build/.github/workflows/ci.yml`), both described at the end.

The two QP solvers implement the direction-finding subproblem of the bundle/gradient-sampling framework. `DualActiveSet` (DAS) is a dual active-set method with an incrementally updated Cholesky factorization and hot-starting; `InteriorPoint` (IPM) is a primal-dual interior-point method with a fraction-to-the-boundary rule and a merit-based step-size sub-QP. Both are mathematically substantial and clearly the result of careful work. Findings below are about robustness, consistency, performance, and maintainability — not the underlying derivations.

---

## Findings that emerged from actually building and testing

Two issues surfaced when I put the code through a fresh CMake build and the test suite on Linux (aarch64, GCC 11, reference netlib BLAS/LAPACK). Both are worth your attention because they're invisible on the author's own machine but real for anyone else.

### A. Unqualified `isnan`/`isinf` — build breaks on strict libstdc++ · **High**
`NonOpt/problems/*.cpp` call `isnan(...)`/`isinf(...)` unqualified (91 call sites across ~14 files, e.g. `MaxQ.cpp:61`, `ChainedCB3_2.cpp:111`, `Test29_22.cpp:99`). Under `<cmath>` with a conforming libstdc++ in strict `-std=c++11` mode, these names live only in `namespace std`, so the compiler rejects them:

```
error: 'isnan' was not declared in this scope; did you mean 'std::isnan'?
```

The core `libNonOpt` library compiled cleanly; only the `problems/` library failed. It builds on macOS/libc++ (where the C names leak into the global namespace), which is why it hasn't bitten you. Fix is mechanical — qualify as `std::isnan` / `std::isinf`. This is precisely the class of defect the proposed CI catches.

### B. Test problem #9 fails in *both* QP solvers on this platform · **Investigate**
Running `testQPSolver` (via `ctest`), 79 of 80 sub-checks pass, but test problem **#9** fails in both solvers:

| Solver | status | iters | KKT error | ‖step‖∞ |
|--------|:------:|:-----:|-----------|---------|
| DAS | 2 (fail) | 240 | **−2.86e+00** | 2.96e+00 |
| IPM | 2 (fail) | 1000 (limit) | 3.67e−02 | 1.00e−01 |

Two things stand out: the DAS **KKT error is negative**, which shouldn't happen for a norm-like optimality residual and points to a sign/label issue or a genuinely corrupted solve; and the IPM simply exhausts its 1000-iteration cap without converging. Because both independent solvers fail on the *same* instance, this is most likely a hard/degenerate QP instance interacting with the (unoptimized reference) BLAS/LAPACK rather than a coincidence.

I can't tell from one platform whether this also fails on your toolchain — which is the point. A cross-platform CI matrix plus a deliberate look at what instance #9 is (and whether the negative DAS KKT value is a reporting artifact or a real failure) would settle it. `tests/testQPSolver.*` is where the instance is defined.

---

## DualActiveSet solver

### 1. `cholesky_tolerance_` is applied inconsistently across the two factorization paths · **Medium (numerical)**
The same option guards two different quantities:

- In `choleskyAugment` (line ~1666) the test is on the **pivot** (squared diagonal): `if (diagonal_squared >= cholesky_tolerance_) factor = sqrt(diagonal_squared); else factor = sqrt(cholesky_tolerance_);`
- In `choleskyFromScratch` (line ~1892) the test is on the **factor diagonal** (already square-rooted by `dpotf2_`): `if (matrix[i*size+i] < cholesky_tolerance_) matrix[i*size+i] = cholesky_tolerance_;`

So one path treats `cholesky_tolerance_` as a floor on `d²` (replacing with `√tol`) and the other as a floor on `d` (replacing with `tol`). For a tolerance like `1e-8`, the effective diagonal floors differ by orders of magnitude (`1e-4` vs `1e-8`) between a from-scratch refactor and an incremental update. Since these two paths are supposed to produce interchangeable factorizations, they should apply the safeguard to the same quantity with the same units. Worth reconciling and documenting which convention is intended.

### 2. Leftover debug string in a production path · **Low, trivial**
`choleskyFromScratch` prints `"ARGH! Replacing Cholesky diagonal of ... with ..."` (line ~1893) to the QP report stream. Recommend a neutral phrasing, e.g. `"Warning: Cholesky diagonal below tolerance; replacing %+.16e with %+.16e."` It's user-visible whenever the safeguard trips.

### 3. `factor_` storage scheme is under-documented and index-fragile · **Medium (maintainability)**
The factor is a flat `double*` indexed `factor_[i * system_solution_length_ + j]`, initially sized `2 * system_solution_length_` and grown by doubling inside `choleskyAugment` via the condition `while ((length - 1) * (system_solution_length_ + 1) >= factor_length_ + 1)`. The original full-size intent (`system_solution_length_ * system_solution_length_`) is preserved only as a commented-out expression at line ~206. The growth condition is correct as far as I can tell but is exactly the kind of hand-derived index arithmetic that's easy to break in a future edit. Two low-risk improvements: (a) add a comment deriving the resize condition, and (b) wrap the factor in a small `struct`/accessor (`R(i,j)`) so the stride logic lives in one place instead of being repeated across ~10 loops. This is the highest maintainability risk in the file.

### 4. Hand-rolled triangular kernels bypass BLAS · **Low (perf, size-dependent)**
`choleskyAugment`/`choleskyDelete` implement the rank-1 update/downdate with explicit scalar loops and Givens-style rotations. That's reasonable for update/downdate (BLAS has no direct primitive), but some inner loops walk columns of a row-major array (`factor_[j * system_solution_length_ + i]`), which is cache-unfriendly, and a temporary `new double[length - index + 1]` is heap-allocated on every augment. For typical bundle sizes this is negligible; if QP dimension grows, a preallocated scratch buffer and attention to access order would help.

### 5. `dpotf2_` (unblocked Cholesky) in the from-scratch path · **Low/Medium**
Confirmed at line ~1888. `dpotf2_` is the level-2 (unblocked) factorization; `dpotrf_` (blocked, level-3 BLAS) is the standard choice and is faster as size grows. Given from-scratch refactors should be the exception, impact is modest, but switching is essentially free.

---

## InteriorPoint solver

### 6. `solveLinearSystem` allocates and copies pivots needlessly · **Low (perf)**
At line ~1143 the routine receives `ipiv[]` but allocates a second `ipiv_local = new int[size]`, factors into it, then copies element-by-element back into `ipiv`. It can factor directly into the caller's `ipiv`, removing an allocation and a loop per solve. Related: the LAPACK workspace query (`lwork = -1`) plus `new double[lwork]` is redone on **every** solve; caching the workspace across the IPM iteration loop avoids repeated allocation.

### 7. Per-iteration `Vector` churn in `computeStepSizes` · **Low (perf)**
`computeStepSizes` constructs ~12 temporary `Vector`s (`r`, `s`, `Qdtheta`, `o`, `p`, `s_plus_p`, `theta1/2`, `v1/2`, `r_dual1/2`) each call, i.e. a dozen heap allocations per IPM iteration. Preallocating these as members (reset each call) would cut allocation traffic on the hot path.

### 8. Unguarded divisions in the step-size sub-QPs · **Medium (numerical robustness)**
The 2-D sub-QP solutions divide by `Q1_22`, `Q2_11`, `Q2_22`, and by Schur complements such as `(Q1_11 - Q1_12*Q1_12/Q1_22)`. When the relevant direction is null (e.g. `du == 0` makes `o` the zero vector and `Q1_22 == 0`), these produce `inf`/`NaN`. IEEE `fmin/fmax` will sometimes mask a `NaN` by returning the other argument, so a bad value can propagate silently rather than failing loudly. Recommend explicit small-denominator guards (or a fallback step) with a clear status return, and possibly an assert in debug builds. This is the IPM analogue of the DAS Cholesky floor.

### 9. `pow(x, 2.0)` for squaring · **Low (perf + accuracy), repo-wide**
`computeStepSizes` (and 30 sites across `src/`) use `pow(..., 2.0)` where `x*x` is both faster and more accurate. E.g. the step-selection compares `pow(r_prim1,2.0) + pow(r_dual1.norm2(),2.0) + ...`. Mechanical replacement with multiplication is a safe, broad win.

### 10. Integer option passed a floating literal · **Low (style)**
`addIntegerOption("QPIPM_iteration_limit", 2e+02, ...)` passes `200.0` where an int is expected (implicit narrowing). Use `200`. (`QPDAS` has analogous integer limits worth checking.)

---

## Cross-cutting: duplication between the two solvers · **Medium (maintainability)**
~24 member functions share identical names and near-identical bodies across the two solvers — `combinationNormInf`, `combinationTranslatedNormInf`, `combinationTranslatedNorm2Squared`, `dualObjectiveQuadraticValue[Scaled]`, `dualSolution[Omega]`, `primalSolution[Feasible][NormInf/Norm2Squared]`, `inexactTerminationCondition`, `checkQuantityCompatibility`, `evaluatePrimalVectors`, `finalizeSolution`, and more. Much of this (the QP data `W/G/b/r`, the combination/primal/dual accessors, the inexact-termination test) is independent of the solve algorithm and could move to the `QPSolver` base class or a shared helper. Consolidating would shrink ~3,700 lines materially and, more importantly, remove the risk of the two copies drifting apart (finding #1 is an example of exactly that kind of drift).

### Minor
- Non-copyable idiom uses private-undefined copy ctor/`operator=`; prefer `= delete` (both solvers).
- Reserved include-guard identifiers `__NONOPTQPSOLVER..._HPP__` (same repo-wide issue as the main review).
- `dualSolutionOmegaLength()` and `vectorListLength()` cast `.size()` to `int` unchecked; fine in practice, but a `static_cast` with an assert documents the assumption.

---

## Suggested QP-specific priorities
1. **Qualify `isnan`/`isinf`** (#A) — unblocks Linux builds; mechanical.
2. **Reconcile `cholesky_tolerance_` usage** (#1) and **guard IPM divisions** (#8) — the two numerical-robustness items.
3. **Investigate test problem #9** (#B) — determine whether it's platform/BLAS sensitivity or a real failure, and whether the negative DAS KKT value is a reporting bug.
4. **Hoist shared code into `QPSolver`** (cross-cutting) — reduces drift risk.
5. Perf polish (#4, #5, #6, #7, #9) as convenient.

---

## Deliverable: CMake build (`nonopt-build/`)
A working CMake build, verified end-to-end on Linux (configure → build → `ctest`), laid out to drop onto the repo root:

```
CMakeLists.txt                      # root: build type, standard, BLAS/LAPACK, install/export
NonOpt/src/CMakeLists.txt           # libNonOpt (static), install + headers
NonOpt/problems/CMakeLists.txt      # libNonOptProblems (static)
NonOpt/exes/CMakeLists.txt          # one exe per .cpp (solveProblem, runExperiment)
NonOpt/tests/CMakeLists.txt         # one test exe per .cpp, registered with CTest
.github/workflows/ci.yml            # CI (see below)
```

Highlights:
- **Defaults to `Release`** (`-O3 -DNDEBUG`) — addresses the biggest performance issue from the first review, where the Makefiles built with `-g` and no optimization.
- `find_package(BLAS/LAPACK REQUIRED)` with a `BLA_VENDOR` hook for OpenBLAS/MKL; portable `${CMAKE_DL_LIBS}` instead of hard-coded `-ldl`.
- `find_package(NonOpt)` support for downstream projects via an exported target (`NonOpt::NonOpt`), plus an `install` target.
- CTest registration of each individual test (their exit codes are meaningful); `testAll` is gated by an output regex until it's changed to `return result;`.

Verified locally:
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release   # Found BLAS / Found LAPACK
cmake --build build -j                            # libNonOpt.a, libNonOptProblems.a, 2 exes, 10 test bins
ctest --test-dir build                            # 8/10 pass (failures = findings A/B above)
```
(During verification I applied a throwaway `std::isnan` fix so the `problems` library would compile; that change is *not* part of the delivered files — it's finding A for you to apply.)

## Deliverable: CI (`.github/workflows/ci.yml`)
- Matrix: Ubuntu (gcc + clang) and macOS (clang); installs BLAS/LAPACK per-OS; configures, builds, runs `ctest`.
- A second job builds via the **original Makefiles** so that path keeps working, running each individual test and scanning `testAll` output for `"failure"`.
- As written, the Linux jobs will fail until finding A is fixed — that's the intended signal.
