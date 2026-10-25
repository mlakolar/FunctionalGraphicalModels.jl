# Milestone 4 FGM migration report

> **Later change (2026-10-25):** FunctionalDataAnalysis.jl now estimates the noise variance with the rotated diagonal over the full interval by default (its issue #17 benchmark). The estimated-noise values below describe this migration as merged; the current value on the test fixture is 0.030864993605722504.

Status: **audit findings F1–F4 addressed; ready for controller validation**.
Base: `origin/bsplinekit`, **6f0eac3d5cbe0fad93e88dcad28f389d89b61495**.
The first follow-up continued from round 1, **f1860c951336be5550b90cf56f7c468bcf0ac682**,
without restarting the migration. No branches, commits, pushes, PRs, or git
configuration changes were made by the developer.
The audit follow-up continued from **d7840ffa8095040868edc1eabc4dca0e5693408a**.

## Audit follow-up (F1–F4)

- **F1:** The migrated neighborhood comparison now reads means, weights,
  eigenvalues and settings from `fit.fpca`, and bases from
  `fit.target_bases`. It no longer fits FDA separately. The base estimator
  does not retain its FPCA fit, so only the base worker repeats its identical
  FPCA call to recover means and eigenvalues; its bases still come from
  `fit.target_bases`. The CSV also compares the actual `fit.coefficients`
  matrices. Each target's basis signs transform both predictor rows and
  response columns. A validation self-test checks independent component
  sign flips, input preservation, and detection of a coefficient change.
- **F2:** The qiao2019 experiment requests `basis=:marginal` and
  `max_k=opts["components"]`, matching the package's fixed truncation.
- **F3:** CSV provenance appends `+dirty` to the recorded HEAD when tracked
  or untracked working-tree changes exist. Source and worker SHA-256 digests
  still identify the estimator source and exact comparison worker. The CSV
  was regenerated after the code and test changes.
- **F4:** `test/doubly.jl` checks dense rejection of each of `bandwidth`,
  `noise_variance`, `grid`, and nondefault `kernel`, plus rejection of mixed
  grid/irregular nodes and incomplete grid data.

This round changes comparison code, experiment FPCA settings and tests. It
does not change package estimator code or expected-value tolerances.

## Author decisions applied (2026-10-24)

1. FDA PR #16 fixed pooled rank upstream. The refreshed FDA `main` is
   **1b75c014778ad374b30fad70fa9083591e1bb546**, source tree
   `f3e6290f513a5c4350ffa7b184eee62a423518f2`. Its pooled `max_k` default is
   `max(1, min(g-2, p*(n-1)-1))`, with explicit cap `min(g, p*(n-1))`.
   FGM calls FDA's pooled route directly. The existing `max_k=ncomponents`
   requests the caller's fixed truncation, including the structural rank cap;
   there is no local rank workaround or replacement eigenstep. A regression
   fixture now retains four pooled components from two nodes with only three
   subjects, reconstructs the observations, and checks the upper cap.
2. FGM adopts FDA's **raw** smoothed covariance diagonal for estimated noise,
   integrated over the middle half, before PSD projection. This resolves the
   previous raw-versus-PSD question. [Qiao, Qian, James and Guo (2020),
   equations (4)–(5) and the sentence after (5)](https://personal.lse.ac.uk/qiaox/qiao.links/2019_DFGMv3.pdf)
   use the local-linear covariance surface and refer noise estimation to
   [Yao, Müller and Wang (2005), Section 2.2, equation (2)](https://utstat.utoronto.ca/fyao/2005-jasa.pdf).
   That estimator and [Li and Hsing (2010), Section 2](https://arxiv.org/pdf/1211.2137)
   estimate noise from unprojected smoothers. The old FGM PSD-diagonal
   subtractor was a package choice. PSD projection increases its diagonal
   and therefore decreases the excess-diagonal noise estimate. FDA retains
   its local-linear surface diagonal; this is not a claim to reproduce Yao
   et al.'s rotated local-quadratic diagonal smoother.

The two-mode fixture in `test/doubly.jl` quantifies both noise changes:

| Noise convention on the same fixture | Estimate | Reason for change |
| --- | ---: | --- |
| Base `6f0eac3`: middle-grid average, PSD subtractor | 0.09556629502099737 | Original FGM implementation |
| Middle-half integral, PSD subtractor | 0.09729427575077350 | Integration instead of a grid-point average |
| Adopted FDA: middle-half integral, raw subtractor | 0.11159325235704108 | Author-approved unprojected-surface convention; **0.0973 → 0.1116** |

See [the diagnostic](../validation/m4_noise.jl) and
[its refreshed CSV](../validation/2026-10-24/m4_noise.csv). The diagonal
smoother and PSD-projected covariance agree between packages to absolute
errors `1.43e-14` and `3.25e-14`, respectively. The change is the noise
formula, not a different covariance smoother. The regression test also checks
the adopted estimate in an independently assembled conditional-score solve.

FDA also rejects requested counts that include modes below its positive-spectrum cutoff,
where old dense FPCA could return numerical null modes. Dense routes retain
FDA's documented `eigen_tol=1e-10`; tuning marks inadmissible ranks as `Inf`
while still evaluating valid lower ranks. No expected-value tolerance was
loosened. No singular zero-noise conditional-score case was encountered.
Such a case remains an error; no ridge or pseudoinverse was introduced.

## Implementation

- FDA owns dense/irregular multivariate data, quadrature, FPCA, prediction,
  reconstruction, PACE smoothers, eigensteps, and conditional scores.
- Graph estimators and node-major `score_covariance` remain in FGM.
  Neighborhood selection still projects **every** node into each target's
  basis. Doubly functional covariance still excludes fitted means and uses
  the uncentered conditional-score second moment with divisor `n`.
- Simulation returns FDA `MultivariateFunctionalData`. Dense estimator and
  smoothing routes validate complete common-grid FDA input.
- The second-difference coefficient P-spline remains in FGM, using BSplineKit.
  Its derivative-penalty counterpart in FDA is a different objective.
  **Candidate FDA addition:** a coefficient-difference P-spline smoother.
  Smoothing settings are retained in each FDA node's metadata, displayed by
  `show(FDA.component(smoothed, j))`, and IDs are retained.
- All graph results retain their graph/solver choices and FDA fits; `show`
  exposes their settings. Tuning retains candidate grids, fold assignments,
  failures, the selected FDA fit, and smoothing metadata.
- Removed exports, without deprecated aliases: `DenseFunctionalData`,
  `SparseFunctionalData`, `FPCAResult`, `fit_fpca`, `transform_fpca`,
  `reconstruct`, `SparseFPCAFit`, `fit_sparse_fpca`. The README lists each.
  The exported `quadrature_weights` is now FDA's function.
- Project sources use sibling FDA and KernelSmoothing `master` checkouts. CI invokes
  `scripts/instantiate.jl`, which resolves KernelSmoothing before FDA in one
  URL-resolution transaction, including on Julia 1.10. ProximalSolvers stays
  transitive. The script restores Project.toml after **all** Pkg operations,
  including resolve and instantiate. No Manifest is part of the diff.
  **FDA and KernelSmoothing must be registered before FGM can be registered.**

## Defaults and explicit choices

| Route | Resolved defaults / required choices | Source and reason |
| --- | --- | --- |
| Dense FPCA | marginal for glasso/neighborhood/doubly; pooled for separable; explicit integer `ncomponents`, `max_k=ncomponents`; no presmoothing | FGM estimator definitions and caller truncation; explicit max_k avoids FDA's automatic smaller cap |
| Dense eigenstep | trapezoidal quadrature; covariance/eigenvalues `n−1`; `eigen_tol=1e-10`, `sign_tol=1e-12`; mean alignment then largest absolute entry | FDA FPCA and multivariate FPCA docstrings, following fdapace 0.6.0 except documented FDA numerical sign/cutoff rules |
| Sparse PACE | required work grid and bandwidth; marginal PACE; Gaussian kernel, optional Epanechnikov; `weighting=:subj`; same bandwidth for mean, covariance, diagonal; `binning=nothing` | Preserve FGM settings, rather than fdapace automatic bandwidth/binning defaults; kernels map to KernelSmoothing objects |
| Sparse spectrum | `max_k=length(grid)`, `eigen_tol=100eps(Float64)`, FDA `sign_tol=1e-12`; all retained positive modes in conditional observation covariance | Preserve FGM full-mode conditional-score algorithm and former relative positive-mode cutoff, with FDA signs |
| Sparse noise | `nothing` estimates; supplied nonnegative scalar or vector passes through per node | Existing FGM interface; adopted FDA middle-half integral using the raw smoothed covariance diagonal, as approved above |
| Score covariance | `center=true`, `corrected=false` (`1/n`) | Existing FGM likelihood scaling; `corrected=true` selects `1/(n−1)` |
| P-spline | required `nbasis`; degree 3, penalty 0, second coefficient differences, trapezoidal loss | Existing FGM cubic unpenalized defaults; FDA derivative penalty does not reproduce the objective |
| Component tuning | `nbasis_grid=5:10`, `ncomponents_grid=1:6`, 5 randomized balanced time folds, supplied RNG or `Random.default_rng()`; degree 3, penalty 0 | Existing FGM implementation of Qiao et al. Section 3.3; exact fold IDs are stored; ties use input order |
| Graphical lasso | required lambda and component count; rho 1, maxiter 10000, abstol 1e-7, reltol 1e-6, threshold 1e-6 | Existing FGM numerical ADMM controls, not paper defaults; objective unchanged |
| Neighborhood | required lambda and dense component count; threshold 0, rule `:or`, center true; maxiter 10000, atol 1e-7, rtol 1e-6, rho 1 | Existing FGM graph/ADMM controls for EJS objective (13) |
| Neighborhood SCV | caller lambda/threshold grids; 5 cyclic subject folds unless explicit fold IDs; OR rule and centered scores | Existing FGM implementation of EJS Algorithm 2; stored fold IDs and candidate convergence |
| Doubly functional | required lambda; ncomponents 2, threshold 1e-6; sparse bandwidth/grid required; dense rejects sparse-only tuning | Existing FGM reconstruction and pointwise CLIME; noise/kernel settings are not silently ignored |
| Separable / joint glasso | required lambda; ncomponents 2, alpha 0.5, threshold 1e-6; joint ADMM rho 1, maxiter 2000, abstol 1e-6, reltol 1e-5 | Existing FGM equal mixture, graph threshold and numerical controls; joint likelihood and penalty unchanged |
| Simulation | n 100, p 50, ncomponents 5, ntime 100, model `:ar2`, noise_sd 0.5, random_diagonal_margin 0.1 | Existing FGM generating-design implementation; settings and caller RNG remain explicit |

The migration adds no new automatic tuning. The additional implementation
choices are explicit max_k overrides, no sparse binning, retention of FGM's
sparse eigenvalue cutoff, per-candidate tuning fits to preserve valid lower
ranks, and metadata storage. There are no outstanding author decisions.
The linked external plan was inaccessible through the browser;
the complete task specification supplied in the request was used.

## Numerical changes and reproducible comparison

1. Dense eigenvalues use `n−1` instead of `n`. The analytic test's old
   `[2, 0.25]` and `[4, 0.09]` are multiplied by `4/3`. Score covariance
   retains divisor `n`. Scores change only by FDA's sign convention, up to
   floating-point eigenvector differences on the fixtures.
2. Sparse estimated noise uses middle-half integration and the approved raw
   covariance subtractor. The table above separates these changes. Increased
   noise changes conditional scores, reconstructed covariances and pointwise
   precisions; these changes are measured below instead of hidden by tolerances.
3. Signs use FDA's mean alignment and largest-entry fallback. Comparisons
   align bases by weighted inner products and transform score-coordinate
   covariance and precision consistently. No tolerance changes hide signs.

Run on one Julia version after instantiation:

```sh
julia +lts --project=. scripts/instantiate.jl
julia +lts --project=. validation/m4_compare.jl results/m4-comparison-new
julia +lts --project=. validation/m4_noise.jl results/m4-noise-new.csv
```

The comparison checks the full base hash, copies the read-only reference into
a temporary directory before resolving/running it, and executes isolated
base/current workers. CSV comments record Julia, dependency tree hashes and revisions, source
and worker digests, HEAD with a dirty-tree marker when applicable, fixture
provenance, and every estimator/FPCA setting. Neighborhood precision is `NA`;
the separate `coefficients_relative_change` CSV column measures its fitted
regression coefficients. Other estimators have `NA` in that column.
Pointwise DFGM edges sum across
time. Relative change is max absolute difference divided by the largest
absolute base entry; scores and matrices are sign-aligned.

All adjacency arrays agree with the base on these fixtures. The sparse rows
exercise both supplied and estimated noise. The direct covariance glasso row
also checks the unchanged low-level solver.

The final column measures regression coefficients for neighborhood selection
and precision matrices for the other estimators.

| Estimator / fixture | Base edges | Migrated edges | Largest relative score change | Covariance change | Precision / coefficient change |
| --- | ---: | ---: | ---: | ---: | ---: |
| Graphical lasso, analytic covariance | 1 | 1 | NA | 0 | 0 |
| Graphical lasso, dense reconstruction fixture | 1 | 1 | 5.23e-16 | 5.74e-16 | 4.83e-14 |
| Neighborhood, fitted target-specific bases and regression coefficients | 2 | 2 | 1.89e-15 | 6.25e-16 | 3.73e-15 |
| Doubly functional, dense | 17 | 17 | 5.23e-16 | 3.46e-15 | 2.78e-15 |
| Doubly functional, sparse, supplied noise | 0 | 0 | 2.53e-15 | 8.55e-15 | 8.52e-15 |
| Doubly functional, sparse, estimated noise | 0 | 0 | 0.00303353 | 0.00507639 | 0.00510229 |
| Partially separable, dense | 1 | 1 | 5.02e-16 | 3.55e-16 | 9.59e-16 |
| Partially separable, pooled rank four with three subjects | 0 | 0 | 1.06e-15 | 5.08e-16 | 7.54e-15 |

Dense eigenvalues agree after the specified `n/(n−1)` adjustment with largest
relative discrepancy at most `2.79e-15`. The machine-readable table is
[validation/m4_comparison.csv](../validation/m4_comparison.csv). The dense
graphical-lasso comparison uses `test/doubly.jl`'s full-rank reconstruction
fixture. `test/fpca.jl`'s two-node fixture shares exactly the same two latent
scores across nodes; round 1 reported that its trial base solver run did not
converge there. No unconverged fit is included in this table.

## Previous follow-up validation

This follow-up used **Julia 1.10.12 only** (`julia +lts`). FDA and
KernelSmoothing were refreshed by `scripts/instantiate.jl`. KernelSmoothing
0.1.0 resolved to tree `3178afc11653d934bdb8aec4a68539fe1473ec46`; dependency
versions, revisions and tree hashes are also recorded in the generated CSVs.
Only the two test files touched in this round were run. No `Pkg.test()` or
docs build was run. Current check outputs:

| Targeted files | Result |
| --- | --- |
| `test/separable.jl` | 65/65 passed, including 10 pooled-rank checks |
| `test/doubly.jl` | 300/300 passed, including estimated-noise and conditional-score checks |
| `validation/m4_compare.jl` | PASS: eight base/current fixture rows, all adjacency arrays equal; both workers completed without warnings |
| `validation/m4_noise.jl` | PASS: both diagnostic rows reproduced; unprojected-surface estimate 0.11159325235704108 |
| `scripts/instantiate.jl` Project.toml contents after execution | byte-for-byte equal to the starting revision, PASS |
| `git diff --check` | PASS |

The targeted test command was:

```sh
julia +lts --startup-file=no --project=. -e 'using FunctionalGraphicalModels, LinearAlgebra, Random, Test; import FunctionalDataAnalysis as FDA; include("test/separable.jl"); include("test/doubly.jl")'
```

Base and current fixture workers are comparisons, not full test runs. No R
runs or golden-value tolerance changes were needed. Round 1's report recorded
passes on Julia 1.13.2 for the other migrated test files; those historical
results do not certify the refreshed dependencies. The controller must run
the complete suite on Julia 1.10 and the current release, including Aqua,
and perform the documentation/CI checks. These acceptance checks remain
pending and are not reported as developer passes.

All three Julia experiment scripts were updated: `qiao2019.jl`,
`qiao2019_comparison.jl`, `other_methods_smoke.jl`. The plotting script uses
CSV only and needs no type migration. `scripts/discovery.sbatch` invokes the
same updated entry points and needs no data changes. It was not submitted.
Historical validation CSVs, manifests, settings and figures were preserved;
they remain evidence about the base implementation, not current validation.
The migration comparison and noise CSVs were refreshed in this round.
Round 1 recorded successful qiao2019 driver `--help` commands and parsing of
the experiment/comparison scripts. This round executed both validation
drivers; the base and current workers emitted no warnings and the
source-digest guard passed.
The repository has Markdown docs and no Documenter project/build entry point;
no claim of a successful documentation build is made.

## Audit follow-up validation

Only **Julia 1.10.12** (`julia +lts`) was used in this audit round. The existing
resolved dependencies were retained; no full suite or docs build was run.

| Check | Result |
| --- | --- |
| `test/doubly.jl` | 306/306 passed, including the six new rejection checks |
| `validation/test_m4_compare.jl` | 6/6 passed for coefficient and projection sign alignment |
| `validation/m4_compare.jl` | PASS: eight rows, all adjacency arrays equal; neighborhood coefficient change 3.734721801467613e-15; both worker logs empty |
| CSV provenance and metrics | PASS: HEAD `d7840ffa8095040868edc1eabc4dca0e5693408a+dirty`; worker/source digests verified; previous metric cells unchanged |
| `git diff --check` | PASS |
| qiao2019 boundary smoke | Five components with six subjects and nine grid points accepted; preprocessing succeeded, two neighborhood fits converged; two glasso fits hit 10000 iterations (exit 1), so this is not a solver-convergence pass |

The targeted test command was:

```sh
julia +lts --startup-file=no --project=. -e 'using FunctionalGraphicalModels, LinearAlgebra, Random, Test; import FunctionalDataAnalysis as FDA; include("test/doubly.jl"); include("validation/test_m4_compare.jl")'
```

The boundary smoke command was:

```sh
julia +lts --startup-file=no --project=. experiments/qiao2019.jl --smoke --n 6 --p 2 --ntime 9 --components 5 --nbasis 7 --nlambda 2 --lambda-min-ratio 0.9 --out results/m4-audit-f2-smoke
```

This small-sample boundary exercises F2: FDA's automatic marginal cap would
be four components. The driver now reaches both graph estimators at the
requested five components. The glasso nonconvergences are retained in the
ignored smoke output, with no preprocessing errors or solver exceptions.
The comparison was regenerated with:

```sh
julia +lts --startup-file=no --project=. validation/m4_compare.jl results/m4-audit-comparison
```

The generated table was copied to `validation/m4_comparison.csv` as the final
tracked artifact update. Its neighborhood settings now show the estimator's
own `max_k=2`. All previously reported metric cells remain unchanged; the new
coefficient column measures the neighborhood output that was missing before.
The full-suite, Aqua, Julia-release, and documentation checks remain the
controller's responsibility.
