# Comparison scope amendment, 25 October 2026

Mladen selected the first 20 fixed replicate IDs per setting for our bias-aware
ACLIME and l2-ACLIME proposals. Default/tuned scaled lasso and every applicable
original-paper comparator retain all 100 replicates per setting. The twelve
settings, seeds, covariance constructions, candidate menus and tuning criteria
are unchanged. The amendment follows observed runtime and completed controls;
partial paper outputs had already been inspected. It is recorded explicitly.

The primary proposal comparison therefore has 240 datasets across twelve
settings. Original methods and scaled lasso have 1,200 intended datasets.
Contrasts involving either proposal use the first 20 paired seeds. Contrasts
between other methods use all 100. Losses are averaged over time before MCSE.

All completed results beyond replicate20 remain retained, including partial
precision checkpoints. They are excluded from the primary proposal summaries
by the fixed replicate rule, rather than chosen according to runtime or loss.
Completed datasets are never refitted for this amendment. A worker that has
already entered scaled-lasso/original-method fitting is allowed to finish.

For incomplete datasets beyond replicate20, the optional driver argument
`--method-subset scaled_original` omits only `aclime` and `l2_aclime`.
Original DFGM-C/DGM-C/Naive-C CLIME comparators remain included. Unchanged,
provenance-matching saved covariances avoid repeating reconstruction. Stopping
and replacing an owned pool requires reconciliation of its terminal state and
completed datasets before any replacement is submitted. Old outputs, sources,
failures and interruption receipts remain intact.

Reporting reads an explicit `comparison-scope.json` at its output root.
Absent that file, legacy reports retain their original scope. Source/protocol
validation recognizes both preserved full-method outputs and the new subset
variant, while checking all rows actually present. Selection/ROC missingness
uses each method's intended replicate set. Missing estimates are not imputed.

This change prioritizes the empirical scaled-lasso comparison. It does not
assert theoretical risk domination or establish Biometrika table replication.
