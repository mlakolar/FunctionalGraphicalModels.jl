# Completed paired implementation comparison

See [the scientific report](../../docs/IMPLEMENTATION_COMPARISON_20261026.md).
All original1200 and new1840 outputs were checked; 240 new paper datasets are
paired to a subset of the original1200, and1600 controls share the prior study.
They must not be added as independent duplicate datasets. Original comparators
use100 seeds per setting; own-method primary contrasts use fixed first20.
Extra72 original proposal datasets remain preserved outside primary summaries.

- original1200: full source/coverage/selection/ROC validation, all selected-fit
  losses/MCSE, source locks, independent replays, raw414 failed-candidate rows.
- returned1200: all default-SZ returned matrices before PD screening, selected
  norm/graph quantities for every applicable comparator, paired counts and PD
  coverage. No baseline fitting was repeated for this supplement.
- spectral1840: strict new-method report with every paired baseline present;
  all failures, norm/graph/PD/ROC/missingness tables and qualified stage timings.
- recipes: exact Slurm sources, intents and returned IDs. Oldv1 returned report
  was cancelled before starting (0seconds); v2 adds all scalar comparator graph
  quantities. Its predecessor remains recorded, not resubmitted uncertainly.

Large CSV/JSON files use .gz. For example:

```python
import csv, gzip
with gzip.open("original1200/validation.json.gz", "rt") as f:
    import json
    checks = json.load(f)
```

The paired figure source reads the original uncompressed remote report. To use
a compressed public CSV with an older plotting script, decompress a copy first.
The current spectral plotting helper also accepts gzip CSVs directly.
All matrix/selection targets are labeled actual/nominal. Undefined AUC or graph
rates stay missing. Methods with no complete selections are visible in coverage
and the returned supplement, even where the old summary has no numeric row.

Original actual/nominal figures and the final paired figure were visually
inspected. The updated spectral paper curve/coverage figure was inspected;
its other three PNGs are byte-identical to previously inspected plots. Both
original full-report replays reuse prior passed unchanged-dataset checks.
See run-receipt.json for actual scope, rows, errors and transfer SHA256.
Legacy reporter receipt wording is retained verbatim; the scientific objective
is the paired implementation comparison stated in the report.
