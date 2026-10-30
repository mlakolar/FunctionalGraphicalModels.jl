# Multiscale numerical regression fixture

`multiscale-id1-full8.f64` is a 50 x 50 x 8 Julia/Fortran-order little-endian
Float64 array (160,000 bytes). It contains the first eight direct covariance
slices from saved simulation dataset 1: constant model, n=200, p=50, m=8,
replicate 1, t=.5. It contains simulated covariance summaries, not subject data.

Source cache: `/project2/mkolar_1314/FunctionalGraphicalModels/20261026-covariance-pilot/outputs/dataset-001/matrices.f64`.
Original 17-slice SHA256: `d05572255e6b4b699e27baef6085d3e82f6ad3c714cc3f093150ed10af1d6b7b`.
Fixture SHA256: `9ccff1f68f46ec818a04ffffe8d3c35d6f073f76f755a0c03578d2b4c01ba71f`.

The test fixes all eight bandwidths, b=.1 and score penalty .05. Baseline
Clarabel 0.11.1 at 1e-7 returns SOLVED candidates that fail the unchanged
unscaled score guard. Fresh standard-tolerance native workspaces reproduce the
issue. Independent JuMP and stricter native solves provide passing comparisons.
The test checks full-column feasibility, objectives against JuMP and retry
receipts, without requiring a particular platform to trigger a numerical failure.

`multiscale-id151-full8.f64` has the same layout for Fourier5, m=8, replicate 1.
Original 17-slice SHA256: `ad65fb72dd35e87f848bdd3b04a0f9dd036fb53f3e5a2992d7659ed6f215c358`.
Fixture SHA256: `5deb4dde268da577d60446d5336ba8a86f62f0b2ad75aa6358b569b337609190`.
At 1e-7, columns 35 and 41 passed the primal guard but underestimated objective
budgets relative to independent JuMP by more than 1e-4 in normalized units.
Native 1e-8 resolves this discrepancy while retaining all acceptance checks.
Both saved fixtures are tested against independent JuMP, with unchanged pilots.
