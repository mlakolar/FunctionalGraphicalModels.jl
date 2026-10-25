"""Validate recorded DFGM coverage and selection consistency without refits.

Run inside Slurm after the producer. Source hashes are checked separately
against the immutable source, not a mutable working tree.
"""
import csv
import json
import math
import pathlib
import sys
import tomllib
from collections import Counter


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def validate(root):
    root = pathlib.Path(root)
    tasks = [line.split("\t")[:2] for line in (root / "tasks.tsv").read_text().splitlines() if line]
    checks = []
    for setting, replicate in tasks:
        folder = root / f"setting-{setting}-rep-{replicate}"
        if not (folder / "complete.toml").exists():
            checks.append(dict(setting=int(setting), replicate=int(replicate), status="incomplete"))
            continue
        settings = tomllib.loads((folder / "settings.toml").read_text())
        receipt = tomllib.loads((folder / "complete.toml").read_text())
        fits = read_csv(folder / "fits.csv")
        counts = Counter(row["method"] for row in fits)
        assert dict(counts) == settings["expected_candidates"], (folder, counts)
        failures = sum(row["status"] != "ok" for row in fits)
        assert failures == receipt["fit_failures"]
        successful = {(row["method"], row["candidate"]): row for row in fits if row["status"] == "ok"}
        metrics = read_csv(folder / "metrics.csv")
        q = len(settings["evaluation_grid"])
        assert len(metrics) == 2 * q * len(successful), folder
        assert all(math.isfinite(float(row[key])) for row in metrics
                   for key in ("operator", "frobenius", "normalized_squared_frobenius")), folder
        with (folder / "selected.csv").open(newline="") as stream:
            selections = {(method, time): candidate for method, time, candidate, status in csv.reader(stream)
                          if status == "selected"}
        assert all((method, candidate) in successful for (method, _), candidate in selections.items())
        for row in read_csv(folder / "auc.csv"):
            method, time = row["method"], row["time_index"]
            if method in ("aclime", "l2_aclime") and len(settings["bias_constant_menu"]) > 1:
                assert int(row["expected_candidates"]) == len(settings["K2_menu"])
                selected = selections.get((method, time))
                if selected is None:
                    assert int(row["successful_candidates"]) == 0
                else:
                    assert float(row["selected_bias_constant"]) == float(successful[method, selected]["bias_constant"])
        assert (folder / "covariances.bin").stat().st_size > 0
        checks.append(dict(setting=int(setting), replicate=int(replicate), status=receipt["status"],
                           candidate_fits=len(fits), fit_failures=failures, metric_rows=len(metrics)))
    (root / "validation.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks))
    return all(row["status"] == "completed" and row["fit_failures"] == 0 for row in checks)


if __name__ == "__main__":
    sys.exit(0 if validate(sys.argv[1]) else 1)
