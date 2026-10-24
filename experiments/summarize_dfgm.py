"""Summarize paired replicate-level results, preserving failures and missingness.

Means and MCSE first average time within a replicate. No time point is treated
as an independent simulation replicate. Python 3.11+; standard library only.
"""
import csv
import math
import pathlib
import statistics
import sys
import tomllib
from collections import defaultdict


def main(root):
    root = pathlib.Path(root)
    groups = defaultdict(list)
    coverage = []
    task_table = root / "tasks.tsv"
    if task_table.exists():
        tasks = [line.split("\t")[:2] for line in task_table.read_text().splitlines() if line]
        folders = [root / f"setting-{setting}-rep-{rep}" for setting, rep in tasks]
    else:
        folders = sorted(root.glob("setting-*-rep-*"))
    for folder in folders:
        receipt = folder / "complete.toml"
        settings = folder / "settings.toml"
        if not receipt.exists() or not settings.exists():
            coverage.append((folder.name, "incomplete", "", ""))
            continue
        s = tomllib.loads(settings.read_text())
        c = tomllib.loads(receipt.read_text())
        coverage.append((folder.name, c["status"], c["fit_failures"], c["elapsed_seconds"]))
        selected = {}
        with (folder / "selected.csv").open(newline="") as f:
            for method, time, candidate, status in csv.reader(f):
                if status == "selected":
                    selected[method, time] = candidate
        values = defaultdict(list)
        with (folder / "metrics.csv").open(newline="") as f:
            for r in csv.DictReader(f):
                if selected.get((r["method"], r["time_index"])) != r["candidate"]:
                    continue
                for loss in ("operator", "frobenius", "normalized_squared_frobenius"):
                    values[r["method"], r["target"], loss].append(float(r[loss]))
        with (folder / "auc.csv").open(newline="") as f:
            for r in csv.DictReader(f):
                values[r["method"], r["target"], "full_auc"].append(float(r["full_auc"]))
        for (method, target, loss), vals in values.items():
            # All intended times must be present and finite for this replicate.
            if len(vals) == len(s["evaluation_grid"]) and all(map(math.isfinite, vals)):
                groups[s["setting_index"], method, target, loss].append(statistics.mean(vals))
    with (root / "summary.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("setting", "method", "target", "quantity", "valid_replicates", "mean", "mcse"))
        for key, vals in sorted(groups.items()):
            se = statistics.stdev(vals) / math.sqrt(len(vals)) if len(vals) > 1 else math.nan
            writer.writerow((*key, len(vals), statistics.mean(vals), se))
    with (root / "coverage.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("task", "status", "fit_failures", "elapsed_seconds"))
        writer.writerows(coverage)
    print(f"Found {len(coverage)} task directories; inspect coverage.csv and per-task fits.csv")


if __name__ == "__main__":
    main(sys.argv[1])
