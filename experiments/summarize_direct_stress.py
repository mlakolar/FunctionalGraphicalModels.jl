"""Validate per-case coverage and summarize held-out-selected stress fits."""
import csv
import math
import pathlib
import statistics
import sys
import tomllib
from collections import defaultdict


def summarize(root):
    root = pathlib.Path(root)
    output, coverage, paired = [], [], []
    for index in range(1, 17):
        folder = root / f"stress-{index}"
        receipt = folder / "complete.toml"
        if not receipt.exists():
            coverage.append((index, "incomplete", 0, 0))
            continue
        completion = tomllib.loads(receipt.read_text())
        settings = tomllib.loads((folder / "settings.toml").read_text())
        with (folder / "stress.csv").open(newline="") as f:
            rows = list(csv.DictReader(f))
        expected = completion["expected_fits"]
        ok = sum(r["status"] == "ok" for r in rows)
        coverage.append((index, completion["status"], len(rows), expected))
        if len(rows) != expected:
            raise ValueError(f"{folder}: {len(rows)} rows, expected {expected}")
        selected = defaultdict(dict)
        for r in rows:
            if r["status"] == "ok" and r["selected"] == "true":
                key = int(r["replicate"])
                if r["method"] in selected[key]:
                    raise ValueError(f"duplicate selected fit: {folder} {key} {r['method']}")
                selected[key][r["method"]] = r
        for method in ("aclime", "l2_aclime", "scaled_lasso"):
            for quantity in ("operator", "normalized_squared_frobenius"):
                vals = [float(methods[method][quantity]) for methods in selected.values() if method in methods]
                if not vals:
                    continue
                if not all(map(math.isfinite, vals)):
                    raise ValueError(f"nonfinite selected loss: {folder}")
                se = statistics.stdev(vals) / math.sqrt(len(vals)) if len(vals) > 1 else math.nan
                output.append((index, settings["n"], settings["p"], settings["degree"], settings["m"],
                    settings["curvature"], settings["noise_variance"], method, quantity,
                    len(vals), statistics.mean(vals), se))
        for rep, methods in selected.items():
            if "aclime" in methods and "l2_aclime" in methods:
                for quantity in ("operator", "normalized_squared_frobenius"):
                    paired.append((index, rep, quantity,
                        float(methods["l2_aclime"][quantity]) - float(methods["aclime"][quantity])))
    for name, header, rows in (
        ("summary.csv", ("configuration", "n", "p", "degree", "visits", "curvature", "noise_variance",
            "method", "quantity", "valid_replicates", "mean", "mcse"), output),
        ("coverage.csv", ("configuration", "status", "recorded_fits", "expected_fits"), coverage),
        ("paired-differences.csv", ("configuration", "replicate", "quantity", "l2_minus_aclime"), paired),
    ):
        with (root / name).open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)
    print(f"Validated {sum(r[1]=='completed' for r in coverage)}/16 completed configurations")


if __name__ == "__main__":
    summarize(sys.argv[1])
