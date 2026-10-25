"""Paired selected-fit contrasts, after averaging time within each replicate.

Run under Slurm. Missing selections never become zeros or unpaired contrasts.
The output gives expected, available and jointly available replicate counts.
"""
import csv
import itertools
import math
import pathlib
import statistics
import sys
import tomllib
from collections import defaultdict


def main(root):
    root = pathlib.Path(root)
    expected = defaultdict(set)
    available = defaultdict(dict)
    applicable = defaultdict(set)
    for line in (root / "tasks.tsv").read_text().splitlines():
        setting, replicate = map(int, line.split("\t")[:2])
        expected[setting].add(replicate)
        folder = root / f"setting-{setting}-rep-{replicate}"
        if not (folder / "complete.toml").exists():
            continue
        settings = tomllib.loads((folder / "settings.toml").read_text())
        applicable[setting].update(settings["expected_candidates"])
        with (folder / "selected.csv").open(newline="") as stream:
            selected = {(m, t): c for m, t, c, status in csv.reader(stream) if status == "selected"}
        values = defaultdict(list)
        with (folder / "metrics.csv").open(newline="") as stream:
            for row in csv.DictReader(stream):
                if selected.get((row["method"], row["time_index"])) != row["candidate"]:
                    continue
                for quantity in ("operator", "frobenius", "normalized_squared_frobenius"):
                    values[row["method"], row["target"], quantity].append(float(row[quantity]))
        for (method, target, quantity), observations in values.items():
            if len(observations) == len(settings["evaluation_grid"]) and all(map(math.isfinite, observations)):
                available[setting, method, target, quantity][replicate] = statistics.mean(observations)
    with (root / "paired-comparisons.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("setting", "first", "second", "target", "quantity", "expected_replicates",
                         "first_valid", "second_valid", "paired_replicates", "mean_difference", "mcse"))
        for setting in sorted(expected):
            for first, second in itertools.combinations(sorted(applicable[setting]), 2):
                for target in ("actual", "nominal"):
                    for quantity in ("operator", "frobenius", "normalized_squared_frobenius"):
                        a = available[setting, first, target, quantity]
                        b = available[setting, second, target, quantity]
                        pairs = sorted(a.keys() & b.keys())
                        differences = [a[r] - b[r] for r in pairs]
                        mean = statistics.mean(differences) if differences else math.nan
                        se = statistics.stdev(differences) / math.sqrt(len(differences)) if len(differences) > 1 else math.nan
                        writer.writerow((setting, first, second, target, quantity, len(expected[setting]),
                                         len(a), len(b), len(pairs), mean, se))
    print("Paired contrasts saved; inspect paired counts and coverage before making comparisons")


if __name__ == "__main__":
    main(sys.argv[1])
