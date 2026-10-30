# Native multiscale path validation, October 27, 2026

Targeted source-level validation: Julia 1.10.12 (Apple ARM), Clarabel 0.11.1,
JuMP 1.32.1, HiGHS 1.26.1, BLAS single-threaded. 139 checks pass.
`targeted.jl` loads the actual precision sources into an isolated module with a
minimal unused functional-data type stub. This verifies the precision solvers
and DirectCovariance wrapper without resolving unrelated functional packages.
It does **not** constitute the complete package suite. The same tests are
included in `test/runtests.jl` for the normal package/CI environment.

Recreate the focused validation in a disposable environment using the supplied
Project and Manifest (rename Manifest-julia-1.10.toml to Manifest.toml there):

```sh
FGM_SOURCE=/absolute/path/to/FunctionalGraphicalModels.jl \
OPENBLAS_NUM_THREADS=1 julia --project=/path/to/disposable/environment targeted.jl
```

The deterministic benchmark uses p=12, three covariance problems (each three
widths), six penalty values and three timed repeats after compilation. Each
run checks budgets against native cold solves and the original JuMP path.
`benchmark.csv` records wall seconds and Julia allocated bytes, not peak RSS.
No simulation data, statistical losses or parameter selection are used.

In a complete package environment:

```sh
OPENBLAS_NUM_THREADS=1 julia --project experiments/benchmark_multiscale_path.jl 12 3
```

The native path builds 12 conic workspaces, updates the matrix 24 times and RHS
180 times; independent native solves build 216. These small-problem timings
must not be extrapolated to p=50/100. A production-size timing canary on
Discovery should precede the next experiment. Native solver reuse changes the
performance implementation, not the estimator objective. Infeasibility remains
visible; the new low-level spectral API does not use a zero fallback.

Full-package CI run 40471342413 failed before tests during resolution:
CovarianceSelection master requires ProximalSolvers 0.1, while ProximalSolvers
master is now 0.2. No dependency compatibility bounds were silently broadened.
The separate `Multiscale solver paths` workflow tests these sources with only
the three solver dependencies on LTS and current Julia. It is an additional
focused check, not a replacement for full-package integration validation.
