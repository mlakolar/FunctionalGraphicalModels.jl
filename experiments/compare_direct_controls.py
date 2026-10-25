"""Pair the two proposals with held-out tuned scaled lasso on identical controls.

Run under Slurm. The original default results remain in the input summaries.
Differences are first-minus-second; their MCSE uses paired replicate differences.
"""
import csv
import math
import pathlib
import statistics
import sys


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


control, calibration, output = map(pathlib.Path, sys.argv[1:4])
output.mkdir(parents=True, exist_ok=True)
differences = []
summary = read(control / "summary.csv")
metadata = {int(row["configuration"]): row for row in summary}
scaled_summary = read(calibration / "summary.csv")
for row in scaled_summary:
    if row["method"] == "heldout_selected":
        merged = dict(metadata[int(row["configuration"])])
        merged.update(row)
        merged["method"] = "scaled_lasso_heldout"
        summary.append(merged)
for index in range(1, 17):
    selected = {}
    for row in read(control / f"stress-{index}" / "stress.csv"):
        if row["status"] == "ok" and row["selected"] == "true":
            selected[row["method"], int(row["replicate"])] = row
    for row in read(calibration / f"calibration-{index}" / "fits.csv"):
        if row["status"] == "ok" and row["selected"] == "true":
            selected["scaled_lasso_heldout", int(row["replicate"])] = row
    for first, second in (("l2_aclime", "aclime"),
                          ("scaled_lasso_heldout", "l2_aclime"),
                          ("scaled_lasso_heldout", "aclime")):
        for quantity in ("operator", "normalized_squared_frobenius"):
            reps = sorted(rep for method, rep in selected if method == first)
            assert reps == sorted(rep for method, rep in selected if method == second)
            values = [float(selected[first, rep][quantity]) - float(selected[second, rep][quantity])
                      for rep in reps]
            assert len(values) == int(metadata[index]["valid_replicates"])
            assert all(map(math.isfinite, values))
            differences.append((index, first, second, quantity, len(values), statistics.mean(values),
                                statistics.stdev(values) / math.sqrt(len(values))))
with (output / "summary.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
    writer.writeheader()
    writer.writerows(summary)
with (output / "paired-comparisons.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(("configuration", "first", "second", "quantity", "paired_replicates", "mean_difference", "mcse"))
    writer.writerows(differences)
print("Paired comparison includes all16 settings and all four methods; no truth-based tuning")
