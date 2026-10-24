import csv
import json
import math
import pathlib
import statistics
import sys

root, original = map(pathlib.Path, sys.argv[1:3])
summary, pairs, checks = [], [], []
for index in range(1, 17):
    folder = root / f"calibration-{index}"
    with (folder / "fits.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 18, (index, len(rows))
    default = {int(r["replicate"]): r for r in rows if r["candidate"] == "1" and r["status"] == "ok"}
    selected = {int(r["replicate"]): r for r in rows if r["selected"] == "true"}
    with (original / f"stress-{index}" / "stress.csv").open(newline="") as f:
        old = {int(r["replicate"]): r for r in csv.DictReader(f) if r["method"] == "scaled_lasso" and r["status"] == "ok"}
    assert len(selected) == 3 and len(default) == 3
    discrepancy = max(abs(float(default[r]["operator"]) - float(old[r]["operator"])) /
        (1 + abs(float(old[r]["operator"]))) for r in default)
    assert discrepancy < 1e-7, (index, discrepancy)
    checks.append(dict(configuration=index, default_replay_relative_error=discrepancy,
        failed_candidates=sum(r["status"] == "failed" for r in rows)))
    for label, values in (("default", default), ("heldout_selected", selected)):
        for quantity in ("operator", "normalized_squared_frobenius"):
            vals = [float(r[quantity]) for r in values.values()]
            assert all(map(math.isfinite, vals))
            summary.append((index, label, quantity, len(vals), statistics.mean(vals),
                statistics.stdev(vals) / math.sqrt(len(vals))))
    for rep in range(1, 4):
        pairs.append((index, rep, default[rep]["normalized_squared_frobenius"],
            selected[rep]["normalized_squared_frobenius"], selected[rep]["candidate"], selected[rep]["lambda0"]))
for name, header, rows in (
    ("summary.csv", ("configuration", "method", "quantity", "valid_replicates", "mean", "mcse"), summary),
    ("paired.csv", ("configuration", "replicate", "default_frobenius", "selected_frobenius", "selected_candidate", "selected_lambda0"), pairs),
):
    with (root / name).open("w", newline="") as f:
        writer=csv.writer(f); writer.writerow(header); writer.writerows(rows)
(root / "validation.json").write_text(json.dumps(checks, indent=2)+"\n")
print("Validated 288 calibration fits; default results replayed on all 48 original datasets")
