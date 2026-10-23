# FunctionalGraphicalModels.jl

[![CI](https://github.com/mlakolar/FunctionalGraphicalModels.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/mlakolar/FunctionalGraphicalModels.jl/actions/workflows/CI.yml)

Estimation and inference for graphical models of multivariate functional data.

Development follows a [staged implementation and reproduction plan](docs/DEVELOPMENT_PLAN.md).
Data use `subjects × nodes × grid points`; FPCA scores use
`subjects × nodes × components`.

Resolve the current ecosystem master dependencies with
`julia scripts/instantiate.jl`, then run `julia --project -e 'using Pkg; Pkg.test()'`.
The generated Manifest records the exact dependency revisions. This project
runs numerical validation on USC Discovery under Slurm.
