# Direct DFGM validation and experiment status

Date: 25 October 2026 (Asia/Dubai). This is a development and pilot record.
It is **not** a successful replication of the Biometrika tables.

## Verified checks

- Discovery job **13383332**: Julia 1.13.2; full package tests **4366/4366
  passed**, including Aqua. Allocation: 4 CPUs, 8 GB; elapsed 195 seconds.
- Job **13383375**: focused tests **53/53 passed**, including the corrected
  overlapping-transition generator. A reduced p=10,n=30 common/irregular
  canary completed its two datasets. One ℓ₂-ACLIME candidate returned
  `ALMOST_OPTIMAL` and was rejected. The driver retained its failure and
  returned nonzero. All other candidate fits completed. Allocation: 8 CPUs,
  16 GB; elapsed 94 seconds. Original output remains retained.
- Job **13383412**: focused tests **53/53 passed** with the revised 1e-7
  solver tolerances. The same common-grid canary realization then completed
  **all 28 candidate fits without failure**. Raw covariance, baseline FPCA,
  direct programs, scaled lasso, scalar CLIME/glasso and group glasso ran.
  A one-replicate pilot of all twelve published n,p/schedule settings is
  running. Its final coverage and result validation remain pending.
- Local checks: Julia parse checks on all new Julia source/driver files;
  Python AST checks; shell syntax checks; `git diff --check`. Scientific
  computation is confined to Discovery Slurm allocations.

The focused tests cover stable binomial inclusion probabilities against exact
finite sums, analytic inverse recovery with vanishing bias, recomputed conic
feasibility, interval symmetrization's norm inequalities, explicit raw diagonal
products, boundary/design validation, scaled-lasso diagonal reduction, empty
statistics, the twelve-setting inventory, valid generated covariance matrices,
reproducibility, and exact transition counts.

## Reproduction inputs and outputs

Research workspace registry:
`.slurm/runs/20261025-direct-dfgm.json`.

Discovery root:
`/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct`.
Its immutable `setup-source`, `canary-source`, and `pilot-source` directories
each retain their source and the resolved Julia Manifest. Submission logs and
the original rejected numerical result are retained at that root. Pilot output
is under `outputs/pilot`; it includes its task table, every worker exit,
settings/seeds, target discrepancies, raw paths, candidate metrics, BIC
selections, observed ROC spans and explicit completion receipts.

The user authorized Discovery execution. The slurm-workflow standing allowance
is at most 8 concurrent CPUs, 16 GB, and two hours with at least one allocation
running, across jobs and retries. The first two jobs consumed 289 seconds;
the pilot requests at most 4500 additional seconds. Queue waiting is excluded.
No larger compute allowance has been assumed.

## Outstanding scientific evidence

The pilot fixes smoothing bandwidth and has a coarse precision path and one
replicate. Paper mode renders 1200 datasets and a denser path with whole-curve
bandwidth CV. The exact author code, integration resolution, overlapping-edge
convention and GGL tuning are unverified. Mean/bandwidth choices, Fourier
projection target discrepancy, boundary/synchronous extensions and invalid BIC
candidates remain explicit in [the protocol](../../docs/DIRECT_DFGM.md).

The controlled stress driver is implemented. Its run and result validation
remain pending. No empirical superiority or published-table agreement is
claimed from the software checks or reduced canaries.
