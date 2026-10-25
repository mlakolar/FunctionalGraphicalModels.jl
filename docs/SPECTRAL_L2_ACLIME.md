# Spectral-calibrated multiscale SOCP

This method implements program (P), Sections 3–8 of the October 25 comparison
note. It is a modification of ℓ₂-ACLIME. It is distinct from the eigenvalue-clipped
Sun–Zhang pipeline in `fit_direct(...; method=:scaled_lasso)`.

`spectral_calibration` constructs the exact dyadic grid augmented by the
covariance-diagonal calibration width. `spectral_covariance` reuses matching
occupancy matrices and computes only missing widths. `fit_spectral_direct`
returns a `DirectPrecisionFit` with method `:spectral_l2_aclime`.
`spectral_l2_aclime_path` reuses each column model across an explicit score menu.
Existing precision methods retain their objectives and calibration. A later
CI compatibility repair only renames three internal covariance-moment functions
(`M0/M1/M2` to `moment0/moment1/moment2`). All frozen experiment sources retain
the original identifiers and exact file hashes; no running source is modified.

```julia
cal = spectral_calibration(N, p, m; L=1.0, C_K=1.0,
    kappa_bound=10.0, h_max=0.25, Aconf=1.0)
built = spectral_covariance(data, grid, cal; boundary=:local_linear)
fit = fit_spectral_direct(built.covariance; L=1.0, C_K=1.0,
    kappa_bound=10.0, h_max=0.25, Aconf=1.0,
    Kscore=2.0, noise_variance=0.5)
```

These example constants are the empirical paper protocol below, not certified
class constants. Cached covariance reuse additionally requires
`cache_mean_matches=true`; the caller must first establish identical
observations and the supplied mean.

On the unit domain, set `I=ceil(log2(N*m))+1`, `Jstar=I+2`,
`x=Aconf*log(p*N)+log(Jstar)`, and `ell=(3log(p)+x)/N`. The grid consists of
`h_max*2^(-i)` for i=0,…,I together with
`h_d=min(h_max,(sqrt(ell)/(m*L))^(1/3))`. For every h,
`u=sqrt(ell)*(1+q)+ell*(1+q^2)` and
`v=sqrt(ell)*q+ell*(1+q^2)`, with `q=1/(m*h)`.
Set `d_k=max(A[h_d][k,k],N^(-2))` and `B_j=2*kappa_bound/d_j`.
Minimize t subject to `norm(beta,1)<=t`, `norm(beta,2)<=B_j`,
`abs((A[h]*beta-e_j-b[h])[k])<=lambda[h,k,j]`, and both
`norm(b[h],2)<=C_K*L*h^2*t` and `norm(b[h],2)<=C_K*L*h^2*B_j`.
The nonunit-noise allowance and its units are stated below.

All widths enter one column SOCP. The column has both an ℓ₂ ball with fixed
radius `B_j=2*kappa_bound/d_j` and two bias constraints, with radii `bbar_h*t`
and `bbar_h*B_j`. The implementation uses `g=b/bbar_h` to scale the cones and
scales score rows by a fixed positive path-wide constant. These are algebraic
transformations of the same constraints. No precision pilot, unknown sparsity,
precision one-norm, oracle width choice or covariance eigenvalue clipping enters.
For numerical boundedness the implementation also sets `t<=sqrt(p)*B_j`.
This preserves the optimization problem's minimum and precision solutions.
Indeed, `norm(beta,1)<=sqrt(p)*B_j` follows from the column ball. If a feasible
point has a larger t, replace it by `sqrt(p)*B_j`: this is at least B_j, so
the second bias constraint implies the first one still holds. No sparsity
bound is used in this argument.

The original interval symmetrization preserves column diagonals. Its output is
radially capped in induced one-norm at the number of fitting subjects. There is
no PSD precision repair or change in the adjacency threshold (1e-6).

Every solved column records its termination, primal/dual status and residuals
for the one-norm budget, column ball, score inequalities and **both** bias bounds.
Only `INFEASIBLE` with a dual infeasibility certificate invokes a labeled zero
column. An inaccurate status, iteration limit or exception is a numerical
failure. It produces no accepted candidate. A zero-column fallback is never
reported as an optimal feasible solution. If beta=b=t=0 satisfies every score
inequality, zero is an analytic optimum with a distinct status; it needs no
numerical solve. No scientific failure is replaced by a baseline graph.

## Units and design extensions

Write T for the observation-domain length. In unit-domain coordinates the
curvature bound is `L_unit=T^2*L`, where L is in physical time coordinates.
Transporting the specified grid gives
`h_d=min(h_max,(T*sqrt(ell)/(m*L))^(1/3))`, `q=T/(m*h)`, and
`bbar_h=C_K*L*h^2`. Here m is the minimum per-coordinate visit count and N counts
subjects, not visits. `C_K` and L are separate inputs; neither is inferred from
the older empirical bias menu. The implementation does not certify L.

The score allowance in original covariance units is
`Kscore*sqrt(d_k*B_j)*(u_h+noise_variance*B_j*v_h)`.
To verify this transport, multiply observations by c. Then A, d, the noise
**variance**, and the operator-bias coefficient multiply by c²; B and beta
divide by c². Both `sqrt(d_k*B_j)` and `noise_variance*B_j` are unchanged, and
`A*beta-e_j-b` is unchanged. This checks the score and constraint units at fixed
widths. It does not assert full scale equivariance of the diagonal floor N⁻²
or the supplied unit-noise grid recipe after changing amplitude units. The
campaign retains its original amplitude units and records every calibration
input. Nonunit-noise use is a protocol extension, not a new theorem here.

The displayed theory assumes independent uniform coordinate schedules, an
interior window, known centering, unit noise and the stated operator-curvature
and asynchronous-envelope conditions. Boundary local-linear weights, shared
stratified schedules, estimated means, and empirical class constants require
separate theory. Their simulation results are labeled empirical extensions.
The proposal has conditional adaptive upper bounds with minimax orders on common
matching parameter ranges. It is not established as unrestricted minimax
optimal, tuning-free, faster or empirically superior.

## Prespecified comparison, before new truth losses are examined

Paper: fixed first 20 replicate IDs in each of the 12 existing settings (240).
Reuse the original seeds, observations, targets and evaluation times. Actual
Fourier-projected precision is primary; nominal dynamic precision is secondary.
Use known-zero centering and the original noise variance 0.5. The exact grid uses
h_max=0.25, Aconf=1, C_K=1, L=1 and kappa_bound=10. L and the condition-number bound
are **empirical sensitivity inputs**, not certified bounds for the paper class.
Kscore=2 is the fixed-input method. It is not certified as sufficiently large by
these simulations. A distinct empirical variant selects over the fixed menu
[0.01, 0.03, 0.1, 0.3, 1, 2] by the existing observed-data pointwise BIC. No L or
condition bound is chosen using truth error. Both variants use the exact grid;
a reduced-grid study, if added, must carry a separate protocol name.

Controls: all 100 recorded replicates in each of 16 configurations (1,600).
Use the same 160 fitting subjects and 40 held-out subjects. The score menu is the
same; selection minimizes held-out covariance–precision residual. The fixed
method remains separate. At the center the precision condition number is at
most 1.6 (Gershgorin), so use the fixed valid bound 2. For curvature parameter c
and the recorded loading norm squared 0.16, use
`L=max(0.01,0.32*c*(1+c/4),0.16*c^2)` and C_K=1. This bounds both pure and mixed
second operator partials over the full unit domain. The stationary case uses
a positive overbound. h_max=0.25 keeps the center window interior. The noise
variances 0.5 and 2 remain recorded extensions. No realized precision condition
number or unknown precision scale is supplied to the estimator.

Report returned-matrix norm errors and PD-screened likelihood selections as
separate populations. Preserve infeasible fallbacks, numerical failures,
undefined graph rates and missing selections. Average times within each seed
before paired differences and MCSE. Report all three matrix losses, graph
recovery, and runtime stages. Per-candidate precision fitting, model reuse,
symmetrization, calibration, covariance construction and scoring are recorded
separately. Menu-wide timing ratios are not algorithm-only speedups.

[Completed new-method fits and dated paired results](SPECTRAL_COMPARISON_20261026.md)
report adverse outcomes and distinguish pending original baselines from
new-method selection failures. Frozen parallel aggregate precision clocks
are qualified; separate controlled benchmarks supply valid wall times.
