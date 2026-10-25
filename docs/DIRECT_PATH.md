# Reusing direct precision programs along a penalty path

`bias_aclime_path` and `fit_direct_path` solve the same convex programs as the
cold `bias_aclime`/`fit_direct` functions. The reference functions are retained.
This is a numerical implementation change, not a new statistical estimator or
an l2 ridge penalty. The l2 method retains its smoothing-bias norm constraint.

```julia
result = fit_direct_path(C, [.03, .1, .3, 1., 2.];
    method=:l2_aclime, pilot_index=2, bias_constant=Lop,
    noise_variance=.5, ntasks=4)
# result.fits[k] is nothing for a failed candidate; inspect result.errors[k].
```

`ntasks` must not exceed available Julia threads. Each time point has private
solver models. Results and diagnostics are assembled in time order. Native
solvers and BLAS use one thread each in the campaign. Reconstruction similarly
partitions independent functional components and uses private CV caches.
The joint GGL objective is unchanged and remains a joint fit across time.

## What is reused [IMPLEMENTED]

For fixed time and bias constant, the first-stage diagonal pilot does not depend
on K2. It is solved once. Each second-stage column model is constructed once;
only normalized constraint right-hand sides change with K2. The subtraction
of e_j is preserved when updating these bounds. The common symmetric interval
minimax-column-l1 program also reuses its model with updated intervals/bounds.

All original objectives, feasibility tests and solver tolerances remain.
Every accepted solve must terminate `OPTIMAL`. An `ALMOST_OPTIMAL` updated
second-stage column gets at most one independent fresh solve with the same
fixed pilot, objective, constraints and tolerances. Successful restarts are
logged in `numerical-restarts.csv`; unresolved errors remain in `fits.csv`.
A failed candidate never supplies precision matrices or selection metrics.
Different numerical paths can select different points in a nonunique optimum;
bitwise equality with independent cold fits is not promised.

The driver records `--path-reuse`, `--ntasks`, source/dependency hashes and
completion outcomes. Completed bias blocks are checkpointed. Covariance caches
must match seed, dimensions, design, grids, noise, Julia/Manifest and the
covariance-relevant source hashes. This permits reuse of expensive reconstruction
from the approved earlier driver without silently reusing an incompatible fit.

## Bounded checks [EMPIRICAL]

Tests cover both norms and LP backends, objective budgets, primal feasibility,
unique-inverse equality, interval symmetrization, input validation and serial
versus threaded assembly. Job 13422768 passed 4,510 package checks. Serial and
threaded FPCA/CV covariances and ranks agree on both common/irregular designs.
On two real saved p=100 time slices, both precision paths agree exactly between
serial and threaded runs at the two tested K2 values.

The cold/path benchmark uses one center covariance slice per dimension, seven
K2 values and fixed bias constant one. On the second measured repeat:

| p | Method | Cold seconds | Reused-path seconds |
|---:|---|---:|---:|
| 50 | ACLIME | 12.26 | 2.06 |
| 50 | l2-ACLIME | 30.07 | 16.60 |
| 100 | ACLIME | 110.62 | 18.55 |
| 100 | l2-ACLIME | 261.80 | 143.50 |

Maximum raw-column and objective-budget differences were below 6e-11 and
7e-12 in these checks. The largest recorded l2 primal residual was about
1.43e-6, within the unchanged budget-scaled 1e-5 acceptance test. These are
bounded implementation comparisons, not campaign-wide speedup estimates or
statistical evidence. [CSV evidence and source receipt](../validation/2026-10-25-campaign/runtime/)
include first-repeat timings and original failure logs.

The full cached pipeline check retains one nonconverged l2 candidate out of
101. The subsequent failure-accounting gate verifies that it is excluded and
counted. Production preserves this policy and reports missing candidates and
invalid positive-definite/BIC selections rather than imputing successful fits.

## FDA default compatibility

The immutable paper campaign retains FDA's earlier surface-diagonal noise
estimator integrated over the middle half of the domain. Its exact covariance
aggregation helpers preserve that convention. Current FDA defaults instead
use a rotated diagonal over the full observed interval. Golden comparisons
therefore request `noise_method=:surface_diagonal` and
`noise_interval=:middle_half` explicitly when those keywords are available.
They retain compatibility with the older locked FDA version and the same
1e-10 agreement tolerance. A change in upstream defaults must not silently
change an active numerical campaign. Current-default noise sensitivity is a
separate source/protocol variant, not a replacement for retained results.
