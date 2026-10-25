# Coverage, selection and independent reporting checks

The paper report uses saved fits. It does not refit models, repair precision
matrices, impute missing selections or append reference ROC endpoints.

## Structural and source checks [IMPLEMENTED]

`validate_dfgm_outputs.py` requires unique, exact keys for every intended fit,
metric, selection and ROC row. Each successful candidate has both target losses
at every evaluation time; failed candidates cannot supply metric rows. The two
targets must have identical BIC values, and selection must choose the first
minimum of the recorded observed-data BIC. All-invalid selections stay missing.

With `--source-locks`, the validator also requires the full 12-by-100 task set,
recognized numerical source hashes and exact Julia/dependency lock. It checks
setting/replicate/seed, n,p,design,schedule, grids, bandwidth/bias/penalty menus,
noise, backend and thread settings against the prespecified campaign protocol.
Source versions are explicit; the optimized shared-grid implementation cannot
silently replace an unrelated dataset or configuration.

ROC areas are reconstructed from the fitted FPR/TPR points. For direct methods,
the BIC-selected bias constant is fixed before varying K2. Duplicate FPR values
use the highest fitted TPR, as in the recorded protocol. Full AUROC requires a
complete candidate path and actual fitted (0,0) and (1,1) endpoints. Partial
area, undefined full AUROC and missing selections remain distinct.

## Saved-array replay [IMPLEMENTED]

`audit_dfgm_saved.jl` reads saved precision paths, covariance inputs and both
precision targets without including the fitting driver or its metric functions.
It computes the operator loss through symmetric eigenvalues, Frobenius losses
through squared entries, support counts through explicit off-diagonal masks,
and Gaussian BIC through a matrix-product trace and Cholesky log determinant.
The campaign's full-matrix threshold count is retained in BIC; this is an explicit
convention, not a switch to counting unique symmetric parameters. The source
reference is Qiao et al. (2020), equation 10 ([DOI](https://doi.org/10.1093/biomet/asz072));
exact author-code equivalence remains unverified.

The validator compares these independently computed rows with the stored report.
The final aggregation checks one completed common-grid and one irregular dataset;
all datasets receive the structural, source, selection and ROC checks. This is a
bounded independent replay, not an independent refit or certification of every
solver optimum.

## Missingness tables [IMPLEMENTED]

- `selection-missingness.csv` separates incomplete datasets, candidate solver
  failures, all-candidate failures and invalid PD/BIC selection times. It counts
  complete, partial and absent time coverage per independent replicate.
- `roc-missingness.csv` separately counts incomplete fitted paths, undefined full
  AUROC and undefined observed area. A complete single-point path still lacks a
  full AUROC.
- `summary.csv` and `paired-comparisons.csv` retain their complete-time replicate
  requirement and report available/paired replicate counts. Those quantities
  describe the available complete cases; differing eligibility must be considered
  before interpreting method comparisons.

## Reduced validation evidence [EMPIRICAL]

Job 13423583 replayed 2,946 metric rows from two completed reduced canaries.
Maximum BIC and loss discrepancies were 1.03e-12 and 2.14e-14. The structural check
rejected fixed-count duplicate metric rows, a nonminimal selected candidate,
invented full AUROC, target-dependent BIC and an unknown source hash. All 244
available production settings at that checkpoint matched the frozen locks.

The reduced irregular canary has 14/21 invalid PD/BIC selections for l2-ACLIME,
5/21 for ACLIME and 12/21 for default scaled lasso, despite zero candidate solver
failures in that dataset. Tuned scaled lasso and the DFGM comparators have no
missing selections there. The common canary has one rejected conic candidate
but complete selected-time coverage. These are coarse three-K2 canary findings,
not results for the full 31-K2 paper study. They remain visible in the evidence.

Job 13423656 checked 17 selection and 34 ROC missingness cells. The saved-array
verifier is byte-identical to the checked version. [Receipts and CSVs](../validation/2026-10-25-campaign/reporting/source-receipt.json)
identify the sources and jobs.

## Final aggregation and retention

Reporter 13423668 replaces the older pending reporter 13423113 and waits for all
300 current producer pools. It uses immutable `paper-report-v3-source` and writes
an attempt-specific directory under `outputs/paper-v4/reports/`. The atomic
`latest-report.json` pointer records its job, directory and exit outcome.
Previous sources, submission records, failures and report attempts are retained.

A successful report job establishes the specified reporting checks. Final
scientific review must still assess failure/missingness patterns, paired
uncertainty, actual versus nominal targets, reference-table discrepancies and
remaining author-code conventions. It does not establish replication merely
because all jobs ended.
