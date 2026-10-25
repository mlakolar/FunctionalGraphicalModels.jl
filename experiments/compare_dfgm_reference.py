"""Compare nominal-target summaries to source-identified author table cells.

No tolerance for declaring replication is implied. Missing/full-AUC values
stay missing; partial ROC areas are never compared to a published full AUROC.
The accompanying fixture JSON identifies the manuscript and version caveat.
"""
import csv
import math
import pathlib
import sys


def main(root, fixture):
    root, fixture = pathlib.Path(root), pathlib.Path(fixture)
    with (root / "summary.csv").open(newline="") as stream:
        summaries = {(r["setting"], r["method"], r["quantity"]): r
                     for r in csv.DictReader(stream) if r["target"] == "nominal"}
    with fixture.open(newline="") as stream:
        reference = list(csv.DictReader(stream))
    with (root / "reference-table-comparison.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("reference_table", "setting", "method", "quantity", "reference_mean",
                         "reference_se", "valid_replicates", "our_mean", "our_mcse",
                         "our_minus_reference", "status"))
        for row in reference:
            value = summaries.get((row["setting"], row["method"], row["quantity"]))
            scale = float(row["display_multiplier"])
            mean, se = float(row["displayed_mean"]) / scale, float(row["displayed_se"]) / scale
            observed = float(value["mean"]) if value else math.nan
            writer.writerow((row["table"], row["setting"], row["method"], row["quantity"], mean, se,
                             value["valid_replicates"] if value else 0, observed,
                             value["mcse"] if value else math.nan, observed - mean,
                             "available" if math.isfinite(observed) else "missing"))
    print("Reference contrasts saved; author-manuscript/final-table equality remains unverified")


if __name__ == "__main__":
    main(*sys.argv[1:3])
