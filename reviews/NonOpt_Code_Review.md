# NonOpt — Code Review & Improvement Suggestions

**Repository:** `frankecurtis/NonOpt` (v2.1) · **Reviewed:** July 12, 2026
**Scope:** code quality & C++, architecture & design, algorithms & numerics, build/docs/tooling
**Basis:** full read of `NonOpt/src/*` plus Makefiles, tests, and examples at the `main` HEAD.

---

## Executive summary

NonOpt is a well-organized, clearly documented package. The strategy-based design is clean and genuinely extensible, the BLAS/LAPACK layer is sensible, and the supporting material (PDF manual, webpage, example problems, AMPL/Python bindings, image-denoising demo) is unusually strong for an academic solver. The suggestions below are mostly about hardening and modernization rather than correctness of the mathematics.

The single highest-value change is trivial: **the library builds with no compiler optimization** (`-g -Wall -std=c++11` only). For a numerical solver this alone can cost a large multiple in runtime. Two other quick wins — a test harness that always reports success, and the absence of CI — are also low-effort, high-impact.

### Priority overview

| # | Finding | Area | Impact | Effort |
|---|---------|------|--------|--------|
| 1 | No optimization flags in any Makefile | Build | **Critical** | Trivial |
| 2 | `testAll.cpp` returns 0 even when tests fail | Tooling | **High** | Trivial |
| 3 | No continuous integration | Tooling | **High** | Low |
| 4 | `Vector::values() const` returns mutable `double*` (const & cache hole) | C++ | **High** | Low |
| 5 | Makefile portability & dependency tracking | Build | Medium | Low |
| 6 | Reserved identifiers (`__NONOPT..._HPP__`) as include guards | C++ | Medium | Low |
| 7 | Exceptions used for normal control flow | Architecture | Medium | Medium |
| 8 | `clock()` for time limits (wrap/multithread issues) | Numerics | Medium | Low |
| 9 | Doc/behavior mismatch: "initialized to zero" constructors don't zero | C++ | Medium | Trivial |
| 10 | Unblocked Cholesky `dpotf2_` on QP matrices | Numerics | Low–Med | Low |
| 11 | Pre-C++11 non-copyable idiom instead of `= delete` | C++ | Low | Low |
| 12 | `normal_distribution` rebuilt each call | Numerics | Low | Trivial |
| 13 | No CMake / install target; duplicate README & LICENSE | Build/Docs | Low | Medium |

---

## 1. Build system

### 1.1 No optimization — **Critical, trivial fix**
Every Makefile (`src/`, `exes/`, `tests/`, `problems/`, `mpc/`) uses:

```make
CXXFLAGS = -g -Wall -std=c++11
```

There is no `-O2`/`-O3` anywhere in the repository (`grep -rn "\-O[0-3s]"` returns nothing), so the shipped library is an unoptimized debug build. For an inner-loop-heavy numerical code this is likely the biggest single performance lever available. Recommended:

- Add a release configuration, e.g. `-O3 -DNDEBUG` (optionally `-march=native` for local builds), and keep `-g -O0` as an opt-in debug target.
- Note that your asserts use the custom `ASSERT_EXCEPTION` macro, so `-DNDEBUG` will *not* silently disable them — you'd need your own guard if you want them compiled out in release. Worth a deliberate decision either way.

### 1.2 Portability and incremental builds — **Medium**
- Compiler is hardcoded (`CXX = g++`). Use `CXX ?= g++` and let `CXXFLAGS`/`LDFLAGS` be overridable so users on clang/macOS/Intel toolchains don't have to edit files.
- `.d` dependency filenames are generated but never used. Add `-MMD -MP` to `CXXFLAGS` and `-include $(depends)` so headers trigger correct recompiles.
- `$(objects): $(sources)` (in `src/Makefile`) makes every object depend on every source, forcing a full rebuild on any single-file edit. Real per-file dependencies (via the `.d` files above) fix this.
- Linker flags are inlined per rule (`-ldl -lblas -llapack`); lifting these into a `LDLIBS` variable makes it easy to swap in OpenBLAS/MKL.

### 1.3 Packaging & adoption — **Low/Medium**
- Consider a CMake build alongside the Makefiles. It would bring `find_package(BLAS/LAPACK)`, an `install` target, `pkg-config`/CMake-config export, Windows/MSVC support, and easy integration for downstream users — a common blocker for adopting academic C++ libraries.
- There is no `install` target and no shared-library option; both help distribution.

---

## 2. Tests & tooling

### 2.1 Test runner masks failures — **High, trivial fix**
`tests/testAll.cpp` accumulates a `result` flag correctly, but `main` ends with:

```cpp
return 0;   // always
```

So the process exit code is `0` even when a suite fails — any CI or shell script checking `$?` will see green regardless. Change to `return result;`. (Also, `testNonOptImplementation()`'s outcome is printed but not folded into `result`.)

### 2.2 No CI — **High**
There is no `.github/workflows`. A minimal GitHub Actions matrix (Ubuntu + macOS, gcc + clang) that installs BLAS/LAPACK, builds, and runs `tests/testAll` would catch regressions and — combined with 2.1 — actually gate on them. This pairs naturally with the `-Wall` you already enable (consider adding `-Wextra`).

### 2.3 Test ergonomics — **Low**
The hand-rolled harness is fine and dependency-free, but exit codes are the only machine-readable signal. Wiring it into CTest (or emitting per-check counts) would make failures easier to triage.

---

## 3. Architecture & C++ quality

The `Strategy` base class with pluggable `DirectionComputation`, `LineSearch`, `QPSolver`, `ApproximateHessianUpdate`, `Termination`, `PointSetUpdate`, and `SymmetricMatrix` implementations is a clean, textbook-appropriate design, and the `Options`/`Reporter`/`Quantities` separation is tidy. Notes:

### 3.1 `Vector::values() const` returns non-const pointer — **High**
In `NonOptVector.hpp`:

```cpp
inline double* values() const { return values_; };
```

A `const Vector&` hands out a writable `double*`. Two problems: it breaks const-correctness, and it lets callers mutate the data **without** going through `valuesModifiable()`, which is the only path that invalidates the cached `norm1_/norm2_/max_/min_/normInf_` values. That's a latent source of stale-cache bugs. Recommended: `const double* values() const` plus a separate non-const overload, and route all writes through the invalidating accessor.

### 3.2 "Initialized to zero" constructors don't zero — **Medium**
`Vector(int length)` and `setLength(int)` are documented as initializing values to zero, but both simply `new double[length]` with no zero-fill. Today `makeNewCopy`/`makeNewLinearCombination` overwrite immediately so nothing breaks, but the contract is false and invites uninitialized reads if anyone relies on it. Either zero-fill (e.g. `dcopy_` from a zero, as the two-arg constructor does) or fix the docs.

### 3.3 No length validation — **Low/Medium**
The default constructor sets `length_ = -1`; `Vector(int)` and `setLength` call `new double[length]` with no check for negative/zero length. A stray call yields UB / `bad_alloc`. A cheap `ASSERT_EXCEPTION(length >= 0, ...)` guards this.

### 3.4 Modernize toward `= delete` and C++17 — **Low**
- Non-copyable classes use the pre-C++11 idiom (private, undefined copy ctor/`operator=`). Since you already require C++11, `= delete` gives clearer compiler diagnostics and states intent (applies across `Vector`, `Strategy`, `RandomNumberGenerator`, etc.).
- Consider raising the baseline to C++17: `std::optional` fits the `Options::valueAs...` lookups, and you'd get `[[nodiscard]]` on status getters, structured bindings, etc.

### 3.5 Manual memory management — **Low**
`grep` finds ~148 `new`/`delete` sites in `src/*.cpp`, mixed with `shared_ptr` ownership. `Vector` in particular manages a raw `double*`. Backing it with `std::vector<double>` (and passing `.data()` to BLAS) would remove the hand-written destructor/`setLength` delete logic and make the class exception-safe with no BLAS downside.

### 3.6 Reserved identifiers — **Medium (portability/standards)**
All 39 include guards use leading double underscores, e.g. `#ifndef __NONOPTVECTOR_HPP__`. Identifiers containing `__` (or a leading underscore followed by a capital) are reserved to the implementation — technically UB and a real clash risk with standard-library internals. Rename to `NONOPT_VECTOR_HPP_`. (Minor related nit: `NonOptVector.hpp`'s closing comment says `__VECTOR_HPP__`, not matching the guard.)

---

## 4. Algorithms & numerics

The mathematical core looks sound: BLAS/LAPACK throughout, `dnrm2` for overflow-safe 2-norms, deterministic seeding for reproducibility, and lazy-cached vector norms are all good choices. Suggestions are about robustness and speed.

### 4.1 Exceptions as normal control flow — **Medium (design)**
The outer loop in `NonOptSolver::optimize` signals *every* exit — including the expected success path — by throwing:

```cpp
THROW_EXCEPTION(NONOPT_SUCCESS_EXCEPTION, "Stationary point found.");
```

caught by a ~20-arm `catch` ladder that maps each exception to a status. It works and is applied consistently, but using exceptions for expected, non-exceptional termination is an anti-pattern: it obscures the control flow and pays throw/unwind cost on the normal path. Consider an explicit status enum returned/broken out of the loop, reserving exceptions for genuine faults. If you keep the pattern, a sentence in the manual explaining the rationale would help future contributors.

### 4.2 `clock()` for time limits — **Medium**
```cpp
if ((clock() - quantities_.startTime()) / (double)CLOCKS_PER_SEC >= quantities_.cpuTimeLimit()) ...
```
`clock()` measures processor time, which (a) can wrap in long runs (`clock_t` is often 32-bit at `CLOCKS_PER_SEC = 1e6`, wrapping near ~36 min), and (b) counts *all* threads' CPU when linked against a multithreaded BLAS, so the "time limit" trips far sooner than wall-clock intent. `std::chrono::steady_clock` gives a monotonic wall-clock limit; expose both if CPU-time is deliberately wanted.

### 4.3 Unblocked Cholesky — **Low/Medium**
`NonOptQPSolverDualActiveSet.cpp:1888` uses `dpotf2_` (the unblocked, level-2 BLAS factorization). For small QP subproblem matrices this is fine, but `dpotrf_` (blocked, level-3) is faster as size grows and is the more standard choice. A size threshold, or just switching to `dpotrf_`, is low-risk.

### 4.4 32-bit BLAS interface only — **Low**
All BLAS/LAPACK prototypes in `NonOptBLASLAPACK.hpp` take `int*`. That precludes ILP64 problems and can silently mismatch a 64-bit-integer BLAS build. Document the LP64 assumption, and consider a single integer typedef to make an ILP64 path possible later.

### 4.5 `normal_distribution` rebuilt per call — **Low**
`RandomNumberGenerator::generateStandardNormal()` constructs a fresh `std::normal_distribution<double>` on every call. Most implementations generate variates in pairs and cache one; rebuilding each time discards that state, roughly doubling engine draws. Make the distribution a data member. (Deterministic seeding to 0 is a nice touch for reproducibility — keep it.)

---

## 5. Documentation & repo hygiene

Strengths worth keeping: the LaTeX/PDF manual, `options.txt`, the project webpage, and the breadth of examples (test problems, AMPL, Python, image denoising) are excellent and above the norm for research solvers.

Small items:
- **Duplicate `README`/`LICENSE`** at the repo root and inside `NonOpt/`. Keep one source of truth to avoid drift.
- **README quick-start.** A copy-paste build recipe with per-OS dependency install (e.g. `apt install libblas-dev liblapack-dev` / `brew install openblas lapack`) and the `NONOPTDIR` setup would lower the barrier for new users, since BLAS/LAPACK aren't bundled.
- A short `CONTRIBUTING` note (build, run `testAll`, style) would help external forks (you already have 9).

---

## Suggested sequencing

1. **Today (minutes):** add `-O3 -DNDEBUG` release flags (1.1); `return result;` in `testAll.cpp` (2.1).
2. **This week (hours):** add a CI workflow (2.2); fix `values() const` const-correctness (3.1); Makefile `CXX ?=` + auto-deps (1.2).
3. **When convenient:** reserved-identifier rename (3.6), `clock()` → `chrono` (4.2), `= delete` modernization (3.4), zero-init doc/behavior fix (3.2).
4. **Larger, optional:** reconsider exception-driven control flow (4.1); add a CMake build (1.3).

*Nothing here indicates a correctness problem in the optimization mathematics — the findings are about performance defaults, safety rails, and modernization.*
