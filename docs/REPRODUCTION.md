# Reproducing the papers

## Implemented estimators

- **FGM:** marginal, quadrature-weighted FPCA; block graphical-lasso objective
  (Qiao et al., equation 9). ADMM solves the same convex objective; it does not
  reproduce the paper's block coordinate-descent implementation or screening.
- **EJS:** target-basis projection, multiresponse block group-lasso, explicit
  threshold and AND/OR rule, lambda-max and selective cross-validation.
  Truncation is supplied; the paper's time-fold basis-dimension CV remains open.
- **DFGM:** local-linear mean and off-diagonal covariance-surface smoothing,
  subject-balanced weights, PSD projection, conditional-expectation scores,
  reconstructed pointwise second moments and exact CLIME. Explicit bandwidth
  and rank; the default noise rule is documented rather than presented as a
  complete fdapace implementation. The `fit_doubly` API takes explicit rank and
  bandwidths. The new `experiments/dfgm2020.jl` driver adds whole-curve bandwidth
  CV, AIC rank selection and pointwise BIC; its candidate grids and remaining
  author-code differences are documented in [the direct DFGM protocol](DIRECT_DFGM.md).
  Gaussian smoothing
  is the default; `kernel=:epanechnikov` satisfies the paper's compact-support
  kernel condition. Gaussian smoothing does not inherit that condition. The
  evaluation grid must cover every observation time. The dense overload
  assumes noise-free or previously smoothed curves.
- **Partial separability:** pooled marginal eigenbasis, score correlations,
  sparse-group joint graphical lasso and union of component supports. Both
  independent-penalty and shared-group limits have analytic tests. Rank and
  penalties are supplied; full FVE/tuning workflows remain to be reproduced.

Dependencies come from `master`: ProximalBase, CovarianceSelection (formerly
CovSel) and its ProximalSolvers dependency (formerly ProximalOPT). B-spline
bases come from the registered BSplineKit.
`FunctionalDataAnalysis` provides local-linear marginal and cross-covariance
surfaces backed by `KernelSmoothing`, which also provides bivariate smoothing.
The DFGM experiment helpers aggregate the exact sufficient statistics of FDA's
raw distinct-visit products and preserve its covariance CV convention. The
new count-normalized direct covariance has a different expectation and is
implemented separately in this package. The
standalone Julia simulation driver adds no data-frame or plotting dependency.

## Seeded Qiao 2019 driver

From the package root, inside a Discovery allocation:

```sh
julia scripts/instantiate.jl
julia --project experiments/qiao2019.jl --smoke --model all --out results/smoke
julia --project experiments/other_methods_smoke.jl results/other-methods
```

`simulate_qiao2019` implements all three Section 5.1 generating structures,
including the original five Fourier coefficients, AR(2) precision weights,
alternating blocks of ten and random edges. The unspecified random-model
positive diagonal and Fourier convention are recorded explicitly.

A bounded paper-size pilot can use:

```sh
julia --project experiments/qiao2019.jl --paper --replicates 1 --p 50 \
  --nlambda 6 --model ar2 --out results/paper-size-pilot
```

`--paper` changes dimensions and default repetition counts. It does **not**
restore the original joint five-fold selection of spline dimension L and
truncation M. This driver fixes L=7 and M (M=5 in paper mode), and compares
FGM with the implemented EJS estimator on the Qiao design. It does not yet
implement ALL/NEVER/HALF/PCA competitors or reproduce the EJS paper's separate
Models A–D. Do not interpret this pilot as the published ROC/AUC table.

The driver writes settings, Julia/source versions, every fitted support,
convergence and KKT diagnostics, and graph metrics. Failures are retained and
produce a nonzero exit. Reference empty/complete graphs are saved separately
and are never substituted for fitted regularization-path endpoints. No AUC
is computed from an incomplete or failed path.

`other_methods_smoke.jl` exercises sparse/dense DFGM and partially separable
fits on a small, explicitly defined finite Gaussian process. Its design is
not the published dynamic-graph or partial-separability simulation. For sparse
DFGM, comparison with the latent pointwise graph is explicitly descriptive:
bounded noisy observations change the reconstruction target.

## Qiao cross-validation and original comparisons

`select_fgm_components` now selects L and M using held-out time points. Each
fold fits B-spline coefficients, means and FPCA from training outcomes only.
The implementation uses a common randomized partition across curves, full-domain
spline support and explicitly documented weighted spline least squares.
Invalid dimensions or rank-deficient candidates are recorded with reasons.

The additional driver implements the five methods in the original Section 5.1:

```sh
julia --project experiments/qiao2019_comparison.jl --p 50 --replicates 1 \
  --nlambda 21 --model ar2 --fgm-rho 0.05 --maxiter 30000 \
  --out results/qiao-comparison
```

FGM uses joint L/M CV over L=5:10, M=1:6. PCA independently chooses L with
M fixed at one. ALL, NEVER and HALF fit scalar glasso on raw noisy observations
at every time point, then respectively take the intersection, union and strict
majority of the estimated supports. Entire diagonals are unpenalized. Penalties
are swept relative to each method's zero-graph bound; the timepoint methods use
one common scalar penalty at each path step. This is an implemented comparison
protocol, not a claim that it reproduces the authors' random-number stream,
unspecified spline conventions or published averaged ROC curves.

Outputs include candidate CV losses and invalidity reasons, selected dimensions,
per-time-point solver diagnostics, metrics and fitted supports. The driver reports
an area over the observed FPR span. It reports a full fitted-path AUC only if
converged fitted empty and complete endpoints are both actually observed.
Reference endpoints remain separate. A single dataset and coarse penalty grid
cannot establish agreement with a 100-replicate published table. The ADMM
parameter `--fgm-rho` changes numerical convergence, not the statistical
objective; difficult small-penalty fits may require a smaller value and a
larger iteration limit. Every such choice is saved in the run settings.

## Discovery execution and receipts

Use `scripts/discovery.sbatch` as a project-specific starting point. Account,
partition, constraint and Julia location were verified for the initial run;
check them again when moving clusters. It requests 2 CPUs, 8 GB and 30 minutes.
It runs tests and small smoke experiments in one allocation with single-threaded
BLAS. The workspace registry `.slurm/runs/20261024-fgm-development.json` records
initial attempts. Exact resolved dependency commits are retained in each
run's generated Manifest. Do not submit full 100-replicate studies without
estimating their cost and respecting the task's approved compute allowance.

## Remaining evidence for full reproduction

1. Complete EJS/partial-separability tuning and competitors. Independently
   validate the Qiao and DFGM comparison conventions against author-generated fixtures.
2. Add the original EJS and partial-separability generating designs. The dynamic
   DFGM design is implemented, but its Fourier-projection target and unspecified
   overlapping-transition conventions still need author-code reconciliation.
3. Compare small numerical fixtures to author implementations where available,
   without copying incompatible licensed code.
4. Run the published dimensions and repetition counts; report Monte Carlo errors.
5. Obtain and verify the real-data inputs and their preprocessing, then reproduce
   the EEG/fMRI examples. No external data are included in this initial package.

These are evidence gaps, not claims about statistical inconsistency of the
implemented estimators. Current solver checks and smoke results are separate
from a successful published-table replication.

An optional scientific plotting helper uses matplotlib without adding Julia
package dependencies. Run it inside a compute allocation with Python 3.11+:

```sh
python3 experiments/plot_qiao.py results/qiao-comparison/metrics.csv results/qiao-roc
```

It exports PNG/SVG figures from the observed fitted points, labels the selected
realization and valid-fit counts, and adds no artificial path endpoints.
