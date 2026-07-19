# CLAUDE.md — NonOpt working brief

This file gives Claude Code the context and task list to improve NonOpt. It
summarizes a prior review; verify each item against the current source before
acting, since line numbers may have shifted.

## Project

NonOpt (Nonconvex, Nonsmooth Optimizer) — C++11 library to minimize locally
Lipschitz (possibly nonconvex/nonsmooth) `f: R^n -> R`. MIT licensed.
Depends on external **BLAS** and **LAPACK** (not bundled).

### Layout
```
NonOpt/src/        core library (Solver, Strategy framework, Vector, QP solvers, ...)
NonOpt/problems/   example/test problems (MaxQ, Chained*, Test29_*, ImageDenoising, ...)
NonOpt/exes/       solveProblem, runExperiment
NonOpt/tests/      testAll + per-component tests (impl lives in the .hpp, thin main in .cpp)
NonOpt/mpc/        image-denoising / QP experiments
NonOpt-Manual/     LaTeX + PDF user manual
docs/              project webpage
```

### Build (existing Makefiles)
```
export NONOPTDIR=$(pwd)
make -C NonOpt/src && make -C NonOpt/problems && make -C NonOpt/tests
cd NonOpt/tests && ./testAll
```
Requires `libblas`/`liblapack` (`apt install libblas-dev liblapack-dev`, or
macOS Accelerate/`brew install openblas lapack`).

## Conventions
- C++11, 2-space indent, `NonOpt` namespace, Doxygen-style comments, brace on
  its own line for functions. Match the surrounding style in any file you edit.
- Vectors/matrices wrap BLAS/LAPACK (Fortran, column-major); NonOpt stores
  row-major and uses the transpose trick — preserve those conventions.
- Non-copyable classes currently use the private-undefined-ctor idiom.

## How to verify changes
After any edit: rebuild `src`, `problems`, and `tests`, then run the individual
test binaries (their exit codes are meaningful) and `testAll`:
```
make -C NonOpt/src && make -C NonOpt/problems && make -C NonOpt/tests
cd NonOpt/tests && for t in testVector testSymmetricMatrix testReporter \
  testQuantities testOptions testPoint testProblem testQPSolver; do ./$t || echo "$t FAILED"; done
./testAll
```
Prefer making one change, rebuilding, and testing before moving on.

---

# Tasks (in priority order)

## 1. Fix unqualified `isnan`/`isinf` — build breaks on strict libstdc++ (HIGH)
`NonOpt/problems/*.cpp` call `isnan(...)`/`isinf(...)` unqualified (~91 sites,
e.g. `MaxQ.cpp`, `ChainedCB3_2.cpp`, `Test29_6.cpp`, `Test29_17.cpp`,
`Test29_22.cpp`). These fail to compile under conforming libstdc++ in strict
`-std=c++11` (Linux/GCC): `'isnan' was not declared in this scope`. They only
compile on macOS/libc++ where the C names leak into the global namespace.
- Replace with `std::isnan` / `std::isinf`; ensure `#include <cmath>` is present.
- Verify: `problems` library compiles under `g++ -std=c++11` on Linux.

## 2. Test suite exit code masks failures (HIGH, trivial)
`NonOpt/tests/testAll.cpp` `main` ends with `return 0;` regardless of the
computed `result`, so CI/scripts always see success. Change to `return result;`
and fold `testNonOptImplementation()`'s outcome into `result` too.

## 3. Enable optimization by default (HIGH, trivial)
All Makefiles use `CXXFLAGS = -g -Wall -std=c++11` — no `-O`. For a numerical
solver this is a large, free speedup.
- Either add an optimized release setting to the Makefiles (`-O3 -DNDEBUG`,
  keep a debug target), or adopt the CMake build (task 9) which defaults to
  Release. Note `ASSERT_EXCEPTION` is a custom macro, so `-DNDEBUG` will NOT
  disable those asserts — decide explicitly whether that's intended.

## 4. Reconcile `cholesky_tolerance_` between the two Cholesky paths (MEDIUM, numerical)
In `NonOptQPSolverDualActiveSet.cpp`:
- `choleskyAugment` floors the **pivot** `diagonal_squared` at `cholesky_tolerance_`
  and sets the factor diagonal to `sqrt(cholesky_tolerance_)`.
- `choleskyFromScratch` floors the **factor diagonal** (already sqrt'd by
  `dpotf2_`) at `cholesky_tolerance_` and replaces with `cholesky_tolerance_`.
These apply the tolerance to different quantities (d² vs d), so the effective
floor differs by orders of magnitude between an incremental update and a
from-scratch refactor. Pick one convention (units), apply it consistently, and
document it. These two factorizations are supposed to be interchangeable.

## 5. Guard divisions in the IPM step-size sub-QPs (MEDIUM, numerical)
In `NonOptQPSolverInteriorPoint.cpp` `computeStepSizes`, the 2-D sub-QP
solutions divide by `Q1_22`, `Q2_11`, `Q2_22`, and Schur complements like
`(Q1_11 - Q1_12*Q1_12/Q1_22)`. A null direction (e.g. `du == 0` -> `o == 0` ->
`Q1_22 == 0`) yields `inf`/`NaN`, and IEEE `fmin/fmax` can silently mask a
`NaN`. Add small-denominator guards with a fallback step and a clear status,
plus a debug assert.

## 6. Investigate QP test problem #9 (INVESTIGATE)
On Linux/aarch64 with reference netlib BLAS/LAPACK, `testQPSolver` problem #9
fails in BOTH solvers: DAS returns status 2 with a **negative** KKT error
(-2.86e+00 — a residual should be >= 0), IPM hits its 1000-iteration limit.
Reproduce on your toolchain/BLAS. Determine (a) whether it fails with your BLAS
or only reference BLAS, and (b) whether the negative DAS KKT value is a
reporting/label bug or a real corrupted solve. The instance is defined in
`NonOpt/tests/testQPSolver.{hpp,cpp}`.

## 7. Reduce duplication between the two QP solvers (MEDIUM, maintainability)
~24 methods are near-identical across `QPSolverDualActiveSet` and
`QPSolverInteriorPoint` (`combinationNormInf`, `combinationTranslatedNormInf`,
`combinationTranslatedNorm2Squared`, `dualObjectiveQuadraticValue[Scaled]`,
`dualSolution[Omega]`, `primalSolution[Feasible][NormInf/Norm2Squared]`,
`inexactTerminationCondition`, `checkQuantityCompatibility`,
`evaluatePrimalVectors`, `finalizeSolution`, ...). Hoist the algorithm-
independent parts (QP data W/G/b/r, combination/primal/dual accessors, inexact
termination) into the `QPSolver` base. This also prevents the kind of drift
behind task 4. Do this carefully, one method at a time, testing between moves.

## 8. Correctness/robustness cleanups in Vector (MEDIUM/LOW)
`NonOptVector.{hpp,cpp}`:
- `double* values() const` returns a mutable pointer from a const object AND
  bypasses the norm-cache invalidation done by `valuesModifiable()`. Provide
  `const double* values() const` + a non-const overload; route writes through
  the invalidating path.
- `Vector(int)` and `setLength(int)` are documented as zero-initializing but do
  not zero the array (`new double[length]`). Either zero-fill or fix the docs.
- No length validation; add `ASSERT_EXCEPTION(length >= 0, ...)`.

## 9. Add the CMake build + CI (MEDIUM)
Prebuilt, verified files exist (from the review) — drop them in and confirm:
```
CMakeLists.txt                     # root: Release default, BLAS/LAPACK, install/export
NonOpt/src/CMakeLists.txt
NonOpt/problems/CMakeLists.txt
NonOpt/exes/CMakeLists.txt
NonOpt/tests/CMakeLists.txt         # per-test CTest registration
.github/workflows/ci.yml            # Ubuntu(gcc,clang)+macOS(clang) CMake+CTest, plus Makefile job
```
Then: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j && ctest --test-dir build --output-on-failure`.
CI (and the CMake build) will fail on Linux until task 1 is done — expected.

## 10. Lower-risk polish (LOW)
- Replace `pow(x, 2.0)` with `x*x` (~30 sites across `src/`, incl. `computeStepSizes`
  and the `Vector` 2-arg constructor) — faster and more accurate.
- `NonOptQPSolverInteriorPoint::solveLinearSystem`: don't allocate `ipiv_local`
  and copy; factor into the caller's `ipiv`. Cache the LAPACK workspace instead
  of re-querying (`lwork=-1`) + `new` on every solve.
- `computeStepSizes`: preallocate the ~12 temporary `Vector`s as members.
- Remove the leftover `"ARGH! Replacing Cholesky diagonal..."` message in
  `choleskyFromScratch`; use a neutral warning string.
- `addIntegerOption("QPIPM_iteration_limit", 2e+02, ...)`: pass `200`, not a
  double literal.
- Consider `dpotrf_` (blocked) instead of `dpotf2_` (unblocked) in
  `choleskyFromScratch`.
- Modernize non-copyable classes to `= delete`.
- Rename reserved include guards `__NONOPT..._HPP__` -> `NONOPT_..._HPP_`
  (leading double-underscore identifiers are reserved to the implementation).
- Timing: `NonOptSolver::optimize` uses `clock()` for the time limit, which can
  wrap and counts all-thread CPU under threaded BLAS; consider
  `std::chrono::steady_clock` for a wall-clock limit.
- Consider whether the exception-driven termination in the main loop
  (`NONOPT_SUCCESS_EXCEPTION`, etc.) should become an explicit status/break;
  at minimum document the rationale.
- De-duplicate the root vs `NonOpt/` `README`/`LICENSE`.

---

## Notes / caveats
- The #9 failure and the isnan build error were observed on Linux/aarch64 with
  unoptimized reference BLAS/LAPACK; confirm behavior on the target toolchain.
- Nothing here indicates an error in the optimization mathematics — items are
  performance defaults, portability, numerical safeguards, and maintainability.
- Suggested commit grouping: (1) isnan + testAll + CMake/CI as a "build/CI"
  commit; (2) cholesky tolerance + IPM guards as a "numerical robustness"
  commit; (3) QP base-class refactor on its own branch; (4) polish separately.
