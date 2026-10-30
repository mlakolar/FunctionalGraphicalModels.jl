# Multiscale ℓ₂-ACLIME paths

## Existing implementations

`bias_aclime_path` / `fit_direct_path(method=:l2_aclime)` solve the two-stage
pilot-based estimator. `spectral_l2_aclime_path` solves the distinct bounded,
spectral-calibrated SOCP. Both retain JuMP models along a score path. HiGHS can
reuse simplex state for the LP stages. Clarabel 0.11's MOI adapter does not
support incremental modification, so modifying a cached JuMP conic model can
recreate the native solver. Retaining the JuMP model alone does not certify
reuse of symbolic factorization or primal/dual iterates.

`multiscale_l2_aclime_path` uses Clarabel's native data updates instead.
A workspace per column retains its sparse structure and symbolic analysis.
Penalty changes update b. Equal-size bandwidth-set changes update A and b.
Explicit zeros preserve the matrix pattern when covariance entries cross zero.
Changing the number of widths or variables builds fresh workspaces. A failed
solve discards that column's workspace before the next attempt. Presolve,
dropping explicit zeros, and chordal decomposition are disabled for updates.
Clarabel restarts its interior-point iterates; this is not iterate warm start.
Pilot LPs are solved once per problem, separately from score penalties.

## Fixed-grid fitting with saved covariance estimates

```julia
using FunctionalGraphicalModels
# C::DirectCovariance contains the saved raw covariance estimates at all widths.
# Indices below refer to C.bandwidths, not bandwidth values.
path = fit_multiscale_l2_path(C, [0.1, 0.2, 0.4, 0.8];
    bandwidth_sets=[[1], [2], [3], [1,2], [2,3], [1,2,3]],
    bias_constants=[0.1, 0.5, 1.0],
    noise_variance=0.5)
# Default pilot is the first width of EACH set; supply pilot_indices to change it.
# Rows identify set × bias constant; columns preserve the supplied penalty order.
path.configurations
path.fits[1,1]       # DirectPrecisionFit or nothing on failure
path.errors[1,1]
path.solver_diagnostics
```

Every width within a set enters the same column program. This is not bandwidth
selection. Each set has its own log(H) confidence resource, matching an
independent `fit_direct` call on that covariance subset. Changing bias constants
changes both the pilot and second-stage envelopes. K1, noise variance and
confidence are also calibration parameters; hold them fixed or explicitly add
separate outer experiment configurations.

The wrapper applies the existing interval-constrained symmetrization and
induced-norm cap n. It preserves the raw covariance input, non-PD successful
estimates, and failed fits. It does not choose parameters or impose a PD screen.
The native solver's epigraph formulation has the same feasible beta/t projection
and objective as the original pilot-based ℓ₂ program. Coefficients need not
coincide for nonunique optima. Status, iterations, setup/solve times and checked
primal residuals accompany each column.

## Explicit spectral calibration or moving covariance problems

```julia
problem = (;matrices=A, u, v, bias, d, kappa_bound=2.0)
r = multiscale_l2_aclime_path([problem], [0.1,0.2,0.4];
    calibration=:spectral, noise_variance=0.5)
columns = r.results[1][1].columns
```

Use `spectral_calibration` / `spectral_covariance` to construct prescribed
spectral resources, widths and diagonal-calibration data. In particular, do not
silently replace the theoretical diagonal width by the first experimental
bandwidth. Each explicit problem may have different matrices, resources, bias
bounds and d. Its kappa bound sets B_j=2*kappa/d_j. Both bias SOCs, the column
ℓ₂ ball and t cap are retained. Input matrices may be indefinite.

The low-level API returns unsymmetrized columns. For a final precision estimate,
apply `interval_symmetrize` and the existing induced-norm cap when required by
the estimator definition. Unlike the older spectral API's explicitly labeled
infeasibility fallback, the new API leaves certified-infeasible columns as NaN
failures. It does not substitute zero. For either calibration, if every width satisfies
`penalty*base[j,h]>=1`, beta=g=t=0 is feasible and globally optimal. This exact
certificate skips the conic solve and reports `termination="ANALYTIC_ZERO"`,
zero iterations and zero setup/solve times. No numerical tolerance is added to
the certificate. The diagnostic `analytic_zeros` counts these column solves.
Skipped points preserve lazy matrix updates even in a nonmonotone penalty path. These failure policies must be distinguished in reports.

The initial native gap/feasibility tolerance is 1e-8, tightened from 1e-7 after
saved-cache objective comparisons exposed inaccurate guard-passing columns.
Pilot LP and symmetrization defaults remain 1e-7. Every original acceptance
constraint is checked at unchanged relative tolerance 1e-5.

A numerically rejected reused solve is rebuilt cold at the same 1e-8 tolerance.
If native recovery fails, an independent JuMP/Clarabel model solves the identical
column program using the original bias vectors b rather than normalized g.
It holds the pilot and every bound/penalty fixed, uses 1e-8, and retries at stricter
1e-9/1e-10 only if needed. This is an equivalent formulation, not covariance or
precision projection. Native SOLVED or JuMP OPTIMAL plus every original guard
is required; ALMOST_SOLVED is never accepted. Certified infeasibility is never
retried or converted into an estimate. Zero bias uses the exact zero auxiliary
witness, which must still pass the original score guard.

`attempt_records` preserves rejected candidates, all statuses, backends, accuracy
settings, residuals, errors and clocks. The column's `backend` identifies the
accepted formulation. `retry_count` counts additional attempts; path diagnostics
include `recovery_attempts`, `cold_restarts`, and `unscaled_solves`. Unknown barrier
iteration counts are explicitly flagged by `iterations_available=false`.

Use `reuse=false` for independent native solves. The returned workspace-build,
matrix-update and RHS-update counts establish which reuse occurred. They do
not establish a speedup; benchmark on the intended dimensions before choosing
cluster allocations.

## Solver validation and next experiment

`test/multiscale_path.jl` compares native reused and independent solves with the
existing JuMP formulations by objective budgets and primal feasibility. It covers
changed covariance sparsity, arbitrary penalty order, changed width counts,
indefinite inputs, infeasibility and recovery, both calibrations, and the final
symmetrization/cap wrapper. Keep solver tests independent of simulation loss.

Before the next statistical campaign:

1. Freeze the covariance cache and target. Save its hash and source revision.
2. Specify singleton bandwidths and multiscale sets as separate estimators.
3. Specify pilot widths, bias envelopes, score paths and noise/condition bounds.
   Separate theory-supplied bounds from empirical fixed constants.
4. Run a small solver-only timing canary with independent-solve checks. Report
   setup, solve, symmetrization, failures and objective/residual differences.
5. Run all fixed configurations on shared replications. Report operator and
   p-normalized squared Frobenius losses, MCSEs, failure rates and paired counts.
   Do not make successful-only averages stand in for unconditional performance.

Solver development lives in FunctionalGraphicalModels.jl because the programs
and their calibration belong to the functional model. CovarianceSelection.jl
continues to provide generic covariance selection methods. No new heavy conic
dependencies are added to that package. No statistical campaign settings or
loss-based parameter selection are embedded in the solver implementation.
