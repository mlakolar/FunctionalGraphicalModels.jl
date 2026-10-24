# FunctionalGraphicalModels.jl

[![CI](https://github.com/mlakolar/FunctionalGraphicalModels.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/mlakolar/FunctionalGraphicalModels.jl/actions/workflows/CI.yml)

Estimation of Gaussian graphical models for multivariate functional data:

- Functional graphical lasso (Qiao, Guo and James, JASA 2019).
- Functional neighborhood selection (Zhao, Zhai, Wang and Kolar, EJS 2024).
- Doubly functional graphical models (Qiao, Qian, James and Guo, Biometrika 2020).
- Partially separable functional graphical models (Zapata, Oh and Petersen,
  Biometrika 2022).

The package uses current `master` versions of **ProximalBase** and
**CovarianceSelection**, with **ProximalSolvers** as a transitive dependency,
and the registered **BSplineKit** for B-spline bases. It supplies weighted FPCA
and sparse functional covariance smoothing locally.

## Setup

Julia 1.10 or newer is required. Resolve ecosystem dependencies explicitly:

```sh
julia scripts/instantiate.jl
```

The generated Manifest records exact dependency commits. The setup script
preserves the package Project file. For an existing resolved environment,
`julia --project -e 'using Pkg; Pkg.instantiate()'` reuses its Manifest.

## Example

```julia
using FunctionalGraphicalModels, Random

sim = simulate_qiao2019(MersenneTwister(24); n=100, p=10)
curves = smooth_bspline(sim.data; nbasis=7)
fgm = functional_glasso(curves, 0.1; ncomponents=5)
ns = functional_neighborhood(curves, 0.1; ncomponents=5, rule=:or)
ps = fit_separable(curves; lambda=0.1, alpha=0.5, ncomponents=5)
dfgm = fit_doubly(curves; lambda=0.1, ncomponents=5)

fgm.diagnostics.converged
all(d -> d.converged, ns.diagnostics)
ps.diagnostics.converged
graph_metrics(fgm.adjacency, sim.adjacency)
```

Data arrays have dimensions `subjects × nodes × grid points`. Score arrays use
`subjects × nodes × components`. Covariance blocks are ordered by node. FPCA
uses trapezoidal quadrature and reuses training means for held-out projection.

`select_fgm_components(sim.data; rng=MersenneTwister(24))` selects spline
dimension and marginal FPCA truncation with held-out-time cross-validation.
The paper-comparison driver includes FGM, ALL, NEVER, HALF and PCA.

The dense DFGM overload expects noise-free or preprocessed curves. For irregular
noisy observations, construct `SparseFunctionalData(times, values)` from matrices
of vectors and call `fit_doubly` with an evaluation grid, bandwidth and component
count. It estimates conditional scores and a graph at each grid point. This
pointwise graph has a different target from a static functional graph.

## Validation and reproduction

All numerical development runs use USC Discovery and Slurm. The suite tests
analytic solutions, projection identities, penalty normalization, conditional
scores, graph metrics and numerical failure handling. See the project-specific
[sbatch script](scripts/discovery.sbatch), [development plan](docs/DEVELOPMENT_PLAN.md)
and [reproduction guide](docs/REPRODUCTION.md). The [validation receipt](docs/VALIDATION.md)
includes a paper-size AR(2) pilot, its ROC figure and adverse findings.

The implementation includes all four estimator families. Published-table and
real-data reproduction remain incomplete. Supplied truncation and tuning choices,
finite ADMM tolerance and functional model assumptions remain explicit. Inspect
convergence diagnostics before interpreting a graph. Sparse CLIME estimates need
not be positive definite; DFGM reports that property at each grid point.
