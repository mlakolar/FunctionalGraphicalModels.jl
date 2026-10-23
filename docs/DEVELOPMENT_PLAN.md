# Development and reproduction plan

The package targets four distinct estimators. Their mathematical targets and
preprocessing must remain explicit. Finite truncation is not automatically an
infinite-dimensional conditional-independence guarantee.

## Stacked pull requests

1. **Functional data foundations.** Common-grid data, quadrature-aware marginal
   and pooled FPCA, out-of-sample projection, reconstruction, B-spline smoothing,
   and covariance with explicit normalization. Use `BSpline.jl` from master.
   Add analytic tests, dependency setup, and this plan.
2. **Functional graphical lasso and neighborhood selection.** Implement Qiao,
   Guo and James's equation (9), with unpenalized diagonal blocks and the ordered
   off-diagonal Frobenius penalty. Use `ProximalBase.LogDetLoss`. Implement the
   EJS function-on-function group-lasso objective, with all predictors projected
   into each target node's basis. Return solver diagnostics and AND/OR graphs.
   Test analytic solutions, KKT conditions, threshold scaling, and invariance.
3. **Doubly functional and partially separable models.** Use current
   `CovarianceSelection.clime` for pointwise DFGM. Add irregular observation
   smoothing and conditional-expectation scores, including measurement noise.
   Implement pooled FPCA and the joint graphical-lasso objective on component
   correlations for partial separability. Check independent limiting cases.
4. **Reproduction tools.** Seeded paper-model generators, graph metrics, regularization
   paths, a small end-to-end reproduction driver, saved settings and diagnostics.
   Run tests and numerical experiments on USC Discovery under Slurm. Separate
   smoke tests from full paper-table replication and document remaining gaps.

Each PR targets the preceding branch; only the first targets `master`. Tests must
pass on the complete stack. Do not change sibling working trees to obtain current
APIs: their master branches contain newer package names and interfaces.

## Paper-to-method map

| Source | Implemented target | Reproduction anchor |
|---|---|---|
| Qiao, Guo, James, *Functional Graphical Models*, JASA (2019), author manuscript supplied by user | Marginal FPCA and block graphical lasso, equation (9) | Section 5.1: five Fourier functions, n=100, p=50/100/150, T=100, noise SD=.5, AR(2), alternating blocks, random graph; 100 replications |
| Zhao, Zhai, Wang, Kolar, EJS 18 (2024), 1042–1129, DOI 10.1214/24-EJS2219 | Target-basis regression and block group-lasso neighborhood selection | Equation (13), Algorithm 1; projection choice, truncation bias and SCV matter |
| Qiao, Qian, James, Guo, Biometrika 107 (2020), 415–431, DOI 10.1093/biomet/asz072 | Conditional reconstruction, pointwise covariance and CLIME | Sparse and dense observation regimes; evolving graph; n=200, p=50/100, 100 replications |
| Zapata, Oh, Petersen, Biometrika 109 (2022), 665–681, DOI 10.1093/biomet/asab046 | Pooled eigenbasis and sparse group joint graphical lasso on correlations | Equation (10); union of component graphs; p=50/100/150, T=30, 100 replications |

Primary sources:

- https://personal.lse.ac.uk/qiaox/qiao.links/research_papers/2017_FGM.pdf
- https://arxiv.org/abs/2105.02487 (published version and author revision)
- https://personal.lse.ac.uk/qiaox/qiao.links/2019_DFGMv3.pdf
- https://arxiv.org/abs/1910.03134

## Design contracts

Dense values have dimensions `n × p × T`; score arrays have dimensions `n × p × M`.
A flattened covariance uses node-major component blocks. FPCA inner products use
trapezoidal quadrature on the supplied, possibly nonuniform grid. All estimated
means come from training data. Solvers expose convergence; a finite iterate is
not silently declared a successful fit. Graph thresholds are explicit.

Use existing packages where their current APIs implement the required primitive:
`BSpline`, `ProximalBase`, `CovarianceSelection` (formerly `CovSel`), and
`KernelSmoothing` when its univariate kernels/smoothers are applicable. Implement
small functional-data components locally where no matching primitive exists.
Do not replace the DFGM covariance surface smoother with independent curve
smoothing or replace the separable joint penalty with independent glasso fits.

## Validation and scope

Run all Julia execution on Discovery. The initial allocation is bounded by the
slurm-workflow standing allowance: aggregate at most 8 CPUs, 16 GB and two hours
of running time. Start with a single small CPU allocation, retain final logs and a
recreation recipe, and record jobs in the workspace's `.slurm/runs` registry.

Initial delivery can establish implemented estimators and validated simulation
infrastructure. Exact published tables require matching all tuning, competitors,
replications and source data. Every experiment must state deviations. Sparse
DFGM reconstruction with bounded measurements targets the reconstructed process;
it need not recover the original pointwise precision. Partial-separability graph
interpretation requires that assumption. No reproduction is certified from a
smoke test.
