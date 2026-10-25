# Exact common-grid mean and covariance calculations

The optional `--common-fast 1` experiment path preserves the local-linear
objectives, fixed-mean whole-curve CV, AIC selection and noise convention of
the frozen paper campaign. It changes their evaluation through sufficient
statistics. The retained `dfgm_common.jl` helper remains the reference.
Irregular observation designs continue to use their existing code.

## Mean curve CV [PROVED: algebraic identity]

Suppose all n subjects have the same T observation times. Let Y contain their
outcomes as T-by-n columns and let H_h be the T-by-T local-linear smoother
matrix at those times. Under the campaign's subject weights, each pooled
observation has the same weight 1/T. A common multiplicative weight cancels
from weighted least squares. Consequently, the fitted values after leaving
subject i out are

```
(n * H_h * mean(Y; dims=2) - H_h * Y[:, i]) / (n - 1).
```

This equals fitting the pooled observations from the other n-1 subjects. The
held-out squared-error average is therefore unchanged. Candidate order and
first-minimum tie selection are preserved. After choosing h, the code calls
FDA's public fixed-bandwidth mean fit on the original observations, preserving
its final values and version-specific metadata. CV scores and assignments
remain recorded. Floating-point arithmetic can perturb near ties; general
bitwise identity of CV scores or selections is not claimed. The numerical
gate required identical bandwidths and ranks on its fixtures. Only complete shared grids and explicit leave-one-curve-out
bandwidth menus are supported by this shortcut.

## Covariance products [PROVED: algebraic identity]

For subject residuals r and Gaussian product weights at query (s,t), define

```
M0(s) = sum_d K((time[d]-s)/h) * r[d]
M1(s) = sum_d K((time[d]-s)/h) * ((time[d]-s)/h) * r[d].
```

Each normal or response moment over distinct visits follows
`sum_{d != e} a[d]*b[e] = sum(a)*sum(b) - sum(a .* b)`.
Let (g0,g1,g2) be the first row of the inverse three-by-three local-linear
normal matrix computed using those distinct-visit moments. The fitted
subject product is

```
g0*M0(s)*M0(t) + g1*M1(s)*M0(t) + g2*M0(s)*M1(t)
- sum_d K_s[d]*K_t[d] * (g0 + g1*u_s[d] + g2*u_t[d]) * r[d]^2.
```

The subtraction removes the same-visit products exactly. Averaging these
fits, or applying the same leave-one-curve deletion formula as above, gives
the original covariance and covariance-CV predictions. The mean is held fixed
through covariance CV, as in the reference FDA path. No low-rank kernel
approximation, binning or altered bandwidth menu is used.

Geometry depends only on observation times, queries and h. Private per-thread
caches reuse it across functional components. At the T-by-T CV queries, the
reference dense pair smoother stores order T^4 entries. The new diagonal
correction stores order T^3 entries and avoids multiplying a T^2-by-T^2
smoother by every subject's product vector. Final covariance surfaces use
the same linear identity applied to `R'R/n`.

The noise estimator remains surface-diagonal integration over the middle
half. Current FDA defaults differ; see [the compatibility note](DIRECT_PATH.md#fda-default-compatibility).
The precision programs, bias/K2 menus, BIC and ROC policies are unchanged.
The alternate implementation and source hashes are explicit in each dataset.

## Bounded numerical evidence [EMPIRICAL]

Discovery job 13423276 compared reference and fast calculations at T=6,10,25,50.
On a regenerated n=200,p=100,T=50 dataset matching the saved frozen input:

- Mean CV took 27.23 s for one reference component and 0.0612 s for the shortcut.
  Final mean values were identical; the largest CV-score difference was1.33e-15.
- Covariance CV/surface fitting took 1.13 s versus 0.155 s for that component.
  The largest surface difference was1.89e-14. Noise, selected bandwidth and
  the tested AIC rank/reconstruction agreed.
- Complete four-thread reconstruction of all 100 components took 10.34 s.
  Its covariance differed from the saved frozen reconstruction by at most
  2.80e-14, and all 100 selected ranks matched.

These are representative checks, not a universal speedup claim. Original
reference timings, source hashes and the retained covariance-only first attempt
are identified in [the evidence receipt](../validation/2026-10-25-campaign/common-fast/source-receipt.json).
The end-to-end common/irregular canary13423280 counted all 101 intended fits:
100 accepted and the same one excluded conic failure. Selection, finite losses,
source correspondence and ROC consistency were checked. This does not certify
zero numerical failures in production.

## Campaign dispatch

After validation, 125 still-pending common-grid pools (500 datasets) moved to
immutable `common-fast-v2-source`. The old pending allocations were cancelled
with a state filter; running jobs retained their original source. Replacement
pools retain four 4-CPU workers and 32 GB, with conservative 3/4/6/6/8-hour limits
for settings 1/2/3/4/5. The original largest common setting was already running.
Seeds, statistical settings and the dependency lock are unchanged.

The journal preserves old/new IDs and canonical output mappings. Dataset
links in `outputs/paper-v4` resolve to the actual implementation attempt.
Aggregation13423113 waits for both original and replacement jobs. Historical
sources, cancelled IDs, failures and partial outputs remain retained.
