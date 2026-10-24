# Direct doubly functional graphical models

This implementation adds direct covariance estimation and three precision
estimators. It preserves `fit_doubly`, the existing Qiao et al. reconstruction
method, as a separate estimator. The new entry point is `fit_direct`.

## Estimators and their scope

**IMPLEMENTED:** bias-aware two-stage ACLIME, multiscale ℓ₂-ACLIME, and scaled
lasso on a projected smoothed covariance. **CONDITIONAL:** the local research
report's theoretical guarantees require its sampling, smoothness, spectral,
pilot, calibration and regular-range assumptions. A numerical solver result
does not verify those assumptions. **EMPIRICAL:** boundary equivalent kernels,
synchronous/stratified schedules, estimated centering, and scaled lasso applied
to this covariance require separate analysis. No optimality or new scaled-lasso
rate is asserted here.

### Raw covariance

Inputs use `FunctionalDataAnalysis.MultivariateFunctionalData`. With independent
uniform observation times on a declared domain, let N_ij count visits in the
target's bandwidth window and π_q be the probability of at least q visits.
The routine uses each curve's own total visit count in this exact binomial
probability. At an interior point on [0,1], define R_ijr=2K((T_ijr-t)/h)Y_ijr,
where K is KernelSmoothing's Epanechnikov kernel. For centered observations,

```
Q_ij = 1{N_ij ≥ 1} sum_r R_ijr / (π_1 N_ij)
H_i,jk = Q_ij Q_ik                                  (j ≠ k)
H_i,jj = 1{N_ij ≥ 2} sum_{r ≠ s} R_ijr R_ijs /
         (π_2 N_ij (N_ij-1))
A_h = sum_i H_i / n
```

Empty windows contribute zero. The diagonal excludes repeated observations;
independent measurement-noise variance therefore needs no diagonal subtraction.
Noise variance still enters the precision estimator's stochastic allowances.
Covariances can be indefinite and are retained as estimated for both ACLIME
methods. Inclusion probabilities use a stable recurrence for the two-visit tail.

The default mean is the explicitly known zero mean. For unknown means, fit on
independent subjects and supply `mean=(j,t)->...`; record the split and fitting
choices. Same-sample mean fitting does not inherit the known-mean theorem.

For a shared random grid, `design=:synchronous` uses distinct-visit pairs for
all entries. Products of marginal window averages would incorrectly reuse
the same visit across nodes. This design has different concentration behavior.
The paper's stratified common grid is an empirical application of this option;
its inclusion probabilities are not the iid binomial probabilities.

`boundary=:error` is the default. `:local_linear` uses analytic truncated
Epanechnikov moments to enforce zero first moment. It can have signed weights.
The existing interior nonnegative-kernel theorem is not asserted at boundaries.

FDA already owns local-linear marginal and cross-covariance smoothing, backed
by KernelSmoothing. It is used for the original reconstruction comparator.
The count-normalized estimator above has a different expectation and proof
structure, so it is implemented here rather than silently replacing FDA's
general-purpose smoother.

### Bias-aware ACLIME

For each column j, first solve an LP minimizing t with -z≤β≤z, z≥0,
t=sum(z), β_j≥0, and

```
|(A_h β - e_j)_k| ≤
  K1 * (max(d_k,d_j)*u_h + noise_variance*v_h) * β_j + b_h*t,
d_k = max(A_h,kk, n^(-2)).
```

Its jth coordinate θ_j supplies the precision-diagonal pilot. The second-stage
LP minimizes the same budget under

```
|(A_h β - e_j)_k| ≤ λ_h,k,j + b_h*t,
λ_h,k,j = K2*sqrt(d_k*θ_j)*(u_h + noise_variance*θ_j*v_h).
```

For normalized domain length one, q_h=1/(m_min*h),
u_h=sqrt(ell)*(1+q_h)+ell*(1+q_h²), and
v_h=sqrt(ell)*q_h+ell*(1+q_h²).
Here ell=(3log(max(p,2))+confidence+log(H))/n. Other domain lengths are
included in q_h. The user supplies b_h=`bias_constant*h²`. The minimum visit
count is a conservative resource input when counts vary. The pilot bandwidth
is explicit. This extends the ACLIME structure; calling the ordinary iid
`CovarianceSelection.aclime` would omit both visit resources and smoothing bias.

### ℓ₂-ACLIME

This is **not an ℓ₂/ridge penalty on the coefficients**. For every supplied
bandwidth h, introduce its own vector b^(h), and solve one SOCP per column:

```
minimize t
|(A_h β - e_j - b^(h))_k| ≤ λ_h,k,j       for all k,h
||b^(h)||₂ ≤ bias_constant*h²*t           for all h
-z ≤ β ≤ z, z ≥ 0, t = sum(z).
```

The same pilot is reused across bandwidths. Every bandwidth is constrained
simultaneously; there is no truth-based bandwidth selection. A finite menu
is supported. To study the report's full adaptive construction, supply its
dyadic grid together with its pilot width and calibrated confidence constants.
The small pilot menu is not presented as the asymptotic adaptive construction.

The ℓ₂ version needs an **operator** bias envelope. An entrywise smoothness
bound alone is insufficient. Both programs have an all-outcome feasible witness
with β=0 and t=max_h 1/b_h. No incorrect universal O(n²) budget is imposed.

### Symmetrization, scaled lasso and solver diagnostics

The corrected common symmetrization solves

```
minimize max_j ||S[:,j]-W[:,j]||₁
S = S',  S_jj = W_jj,
min(W_ij,W_ji) ≤ S_ij ≤ max(W_ij,W_ji).
```

It retains pairwise intervals and serves both induced-norm and squared
Frobenius error transfers. A final induced-ℓ₁ cap at n matches the report's
bad-event control. Raw W is retained separately. Positive definiteness is
reported, never assumed or silently repaired in the precision matrix.

For scaled lasso, eigenvalues of A_h below the explicit positive `psd_floor`
are raised to that floor. Its operator displacement is saved. The resulting
Gram matrix is passed to `CovarianceSelection.slassoinv` with n **subjects**,
not n times the visit count. Its regression and residual-variance normalization
come from Sun and Zhang (2013). The common interval LP replaces that package's
smaller-entry symmetrization. The paper's iid theory does not automatically
cover a smoothed, projected Gram matrix. Defaults use its quantile penalty;
`lambda0` may also be supplied explicitly.

JuMP/Clarabel solves the LPs and SOCPs. Only `OPTIMAL` solutions with primal
values are accepted. Second-stage constraint violations are recomputed from
the returned primal variables. Failed programs throw; experiment drivers retain
their error records. There is no fallback to a zero or oracle graph.

## Minimal example

```julia
using FunctionalGraphicalModels

# data is centered FDA raw functional data with independent uniform times.
C = direct_covariance(data, [.4,.5,.6]; bandwidths=[.1,.2,.3])
acl = fit_direct(C; method=:aclime, pilot_index=2, bias_constant=L)
l2  = fit_direct(C; method=:l2_aclime, pilot_index=2, bias_constant=Lop)
sl  = fit_direct(C; method=:scaled_lasso, pilot_index=2, bias_constant=Lop)
```

The illustrative grid above is interior only if every width stays in the
declared sampling domain. L and Lop denote actual valid envelopes, including
kernel constants; the API cannot infer their validity from observations.

## Biometrika comparison design

`dfgm2020_settings()` enumerates all twelve settings in Section 4:

| Design | p | Schedule | n | Published replicates |
|---|---|---|---|---|
| common | 50,100 | T=10,25,50 | 200 | 100 per setting |
| irregular | 50,100 | T=6:9,20:30,50 | 200 | 100 per setting |

The generator follows Section 4 and Appendix A: .25 baseline diagonal,
40 initial edges and 10 changes per type for p=50, 160/40 for p=100, strengths
in [.1,.3], .5-long linear transitions, Gaussian temporal correlation with
scale one, noise variance .5, and five orthonormal Fourier functions. For
unpublished smaller canary dimensions, edge counts scale quadratically in p.

Unspecified conventions are explicit: trapezoidal integration at a chosen
resolution, constant/sin1/cos1/sin2/cos2 ordering, and replacing an unfinished
growing ramp if that edge is selected to decay. A currently decaying edge is
not selected again before its decay ends. This preserves the prescribed number
of new changes while resolving overlapping transitions. It is not a verified
copy of an author implementation.

**Target check:** Appendix A projects a nonseparable coregionalization kernel
onto five Fourier functions. The generated process's pointwise covariance is
computed from the resulting coefficient covariance. Its inverse can differ
from the nominal dynamic precision used to construct the kernel. Both targets
and their discrepancy are saved. A dense actual projected precision can make
graph ROC undefined because no null edges remain. This is retained as a finding.

Common-grid comparators are DFGM-C/G, DGM-C/G, GGL and Naive-C/G. Irregular
comparators are DFGM-C/G. DFGM reconstructs marginal FDA PACE scores; AIC selects
rank up to five. Pilot mode fixes its bandwidth at .25. Paper mode supplies
whole-curve CV candidates [.1,.2,.3,.4] for mean and covariance smoothing.
This explicit candidate grid is an implementation convention. Covariance CV
can be expensive because FDA evaluates held-out raw products. DGM uses Gaussian
smoothing with the stated bandwidth order and proportionality constant one;
the paper's exploration of constants up to three is not reproduced by that
single constant. GGL uses the group-only joint graphical-lasso objective and
the package's ADMM solver. These distinctions remain reproduction gaps.

Direct methods use known zero centering, bandwidths [.15,.25,.4], and a common
fixed K2 menu. Their bias envelope is exploratory, not certified from truth.
Scaled lasso uses its fixed default quantile penalty, so it has one fitted
operating point rather than an invented regularization-path AUROC.

Precision tuning uses the paper's pointwise BIC on observed covariance inputs.
Indefinite precision candidates receive infinite BIC. If no positive-definite
candidate exists, selection is missing. Every path loss is also retained so
this selection failure remains visible. No latent target enters selection.

`metrics.csv` stores raw operator and Frobenius norms, normalized squared
Frobenius loss, BIC and graph operating points. Time points are averaged within
a replicate before computing Monte Carlo standard errors across replicates.
Summaries require all intended time points for a replicate. Full path AUROC is
missing unless fitted empty and complete endpoints are observed and every
candidate succeeded. Partial observed areas remain labeled separately. Raw
numbers are not multiplied by ten; Tables 2–3 multiply by ten for formatting.

## Controlled additional simulations

`experiments/direct_stress.jl` uses a sparse star precision of degree 2 or 10
at t=.5, n=200, p=50, m=8 or 30, noise variance .5 or 2, and curvature 0 or 4.
Its process is Z+(1+c(t-.5)²)vU with covariance of Z equal to Σ₀-vv'. Thus the
target is exactly Σ₀ and its operator smoothing bias is exactly
`((1+c*h²/5)^2-1)*norm(v)^2`. This checks bias handling and includes a stationary
control. A subject split (160 training/40 validation) selects the common K2
menu by covariance-precision residual, before evaluating latent target loss.
An exact valid bias envelope is known from this controlled construction and
reported as such; this experiment does not demonstrate adaptation to unknown
smoothness. Additional recommended checks are integration-grid refinement
(101 versus 201), increased Monte Carlo replicates, nonzero estimated means,
and sensitivity to the scaled-lasso covariance PSD displacement.

## Discovery execution and evidence

All scientific computation runs under Slurm. The local task generator only
renders commands. Example from a staged package:

```sh
python3 experiments/make_dfgm_tasks.py /absolute/run/pilot.tsv --mode pilot
sbatch --parsable --job-name=UNIQUE_RUN_TAG \
  scripts/discovery_dfgm.sbatch /absolute/source /absolute/run/pilot.tsv /absolute/run/output
```

`--mode paper` renders 1,200 tasks with 31 penalty candidates. It does not
submit them or imply that their compute cost is approved. Paper mode is a
reproduction attempt with the explicit conventions above, not a certification
of the published tables. The default launcher requests 8 CPUs, 16 GB, 75
minutes; a larger campaign needs measured throughput and an authorized budget.
GNU Parallel launches one-CPU exclusive `srun` workers with single-threaded
BLAS. Source snapshots and resolved Manifests belong to each immutable attempt.

The launcher records every worker exit, preserves individual failures and runs
the inexpensive summarizer even if workers fail. The output task table defines
expected coverage. A new attempt uses a fresh output directory. Existing
completed or partial settings are not silently overwritten. Record exact jobs,
source/dependency receipts, seeds, outputs, checks and consumed running time in
the workspace's lightweight `.slurm/runs/` registry. No cleanup is automatic.

See `validation/2026-10-25-direct/REPORT.md` for checks and results actually
performed. Reduced canaries and one-replicate pilots are not paper replication.

## Sources

- Qiao, Qian, James and Guo (2020), *Doubly functional graphical models in high
  dimensions*, Biometrika 107, 415–431, DOI
  [10.1093/biomet/asz072](https://doi.org/10.1093/biomet/asz072).
  [Author manuscript and supplement](https://personal.lse.ac.uk/qiaox/qiao.links/2019_DFGMv3.pdf),
  Section 2.4, Section 4 and Appendix A. Published-table agreement remains to
  be established; the author manuscript is the implementation source used.
- Sun and Zhang (2013), *Sparse matrix inversion with scaled lasso*, JMLR 14,
  3385–3418, [paper](https://jmlr.org/papers/volume14/sun13a/sun13a.pdf),
  equations (4), (7), (9). The package uses the interval restriction described
  above in addition to the paper's symmetric minimax-column LP.
- Local method derivation: `notes/20261015/results_draft.md`, Sections 2–4;
  `notes/20261019/l2_aclime_proofs.md`; and the corrected feasibility budget
  and common symmetrization in `notes/20261019/l2_aclime_audit.md`, A4/A8,
  in the research workspace. These unpublished notes are not bundled here.
- [Clarabel settings](https://clarabel.org/stable/api_settings/) and
  [JuMP conic constraints](https://jump.dev/JuMP.jl/stable/manual/constraints/).
