"""Expose replicate/time selection missingness separately from solver failures.

Run after validate_dfgm_outputs.py under Slurm. The summary's complete-time
replicate requirement stays unchanged; no estimate or loss is imputed.
"""
import csv
import json
import math
import pathlib
import sys
import tomllib
from collections import Counter, defaultdict
from dfgm_scope import included, load_scope


def main(root):
    root = pathlib.Path(root)
    scope = load_scope(root)
    checks = {(r["setting"], r["replicate"]): r
              for r in json.loads((root / "validation.json").read_text())}
    totals = defaultdict(Counter)
    roc_totals = defaultdict(Counter)
    common = ("aclime", "l2_aclime", "scaled_lasso", "scaled_lasso_tuned", "DFGM_C", "DFGM_G",
              "DGM_C", "DGM_G", "GGL", "Naive_C", "Naive_G")
    for line in (root / "tasks.tsv").read_text().splitlines():
        setting, replicate = map(int, line.split("\t")[:2])
        methods = common if setting <= 6 else common[:6]
        methods = tuple(m for m in methods if included(scope, m, replicate))
        check = checks[setting, replicate]
        for method in methods:
            totals[setting, method]["intended_replicates"] += 1
            for target in ("actual", "nominal"):
                roc_totals[setting, method, target]["intended_replicates"] += 1
        if check["status"] == "incomplete":
            for method in methods:
                totals[setting, method]["incomplete_datasets"] += 1
                for target in ("actual", "nominal"):
                    roc_totals[setting, method, target]["incomplete_datasets"] += 1
            continue
        folder = root / f"setting-{setting}-rep-{replicate}"
        settings = tomllib.loads((folder / "settings.toml").read_text())
        q = len(settings["evaluation_grid"])
        with (folder / "fits.csv").open(newline="") as stream:
            fits = list(csv.DictReader(stream))
        succeeded = Counter(r["method"] for r in fits if r["status"] == "ok")
        failed = Counter(r["method"] for r in fits if r["status"] == "failed")
        for method in methods:
            summary = totals[setting, method]
            missing = check["missing_selected_times"][method]
            summary["completed_datasets"] += 1
            summary["evaluated_times"] += q
            summary["missing_selected_times"] += missing
            summary["failed_candidate_fits"] += failed[method]
            summary["complete_selection_replicates"] += missing == 0
            summary["partial_selection_replicates"] += 0 < missing < q
            summary["no_selection_replicates"] += missing == q
            summary["all_candidates_failed_replicates"] += succeeded[method] == 0
            summary["invalid_PD_or_BIC_times"] += missing if succeeded[method] else 0
        with (folder / "auc.csv").open(newline="") as stream:
            auc = list(csv.DictReader(stream))
        grouped = defaultdict(list)
        for row in auc:
            grouped[row["method"], row["target"]].append(row)
        for method in methods:
            for target in ("actual", "nominal"):
                rows = grouped[method, target]
                assert len(rows) == q
                summary = roc_totals[setting, method, target]
                undefined = sum(not math.isfinite(float(r["full_auc"])) for r in rows)
                summary["completed_datasets"] += 1
                summary["evaluated_times"] += q
                summary["complete_full_auc_replicates"] += undefined == 0
                summary["partial_full_auc_replicates"] += 0 < undefined < q
                summary["no_full_auc_replicates"] += undefined == q
                summary["undefined_full_auc_times"] += undefined
                summary["undefined_observed_area_times"] += sum(not math.isfinite(float(r["observed_area"])) for r in rows)
                summary["incomplete_path_times"] += sum(r["successful_candidates"] != r["expected_candidates"] for r in rows)
                summary["complete_path_undefined_auc_times"] += sum(r["successful_candidates"] == r["expected_candidates"] and
                    not math.isfinite(float(r["full_auc"])) for r in rows)
    fields = ("intended_replicates", "completed_datasets", "incomplete_datasets", "evaluated_times",
              "complete_selection_replicates", "partial_selection_replicates", "no_selection_replicates",
              "all_candidates_failed_replicates", "missing_selected_times", "invalid_PD_or_BIC_times",
              "failed_candidate_fits")
    with (root / "selection-missingness.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("setting", "method", *fields))
        for (setting, method), counts in sorted(totals.items()):
            assert counts["completed_datasets"] + counts["incomplete_datasets"] == counts["intended_replicates"]
            assert (counts["complete_selection_replicates"] + counts["partial_selection_replicates"] +
                    counts["no_selection_replicates"]) == counts["completed_datasets"]
            writer.writerow((setting, method, *(counts[f] for f in fields)))
    fields = ("intended_replicates", "completed_datasets", "incomplete_datasets", "evaluated_times",
              "complete_full_auc_replicates", "partial_full_auc_replicates", "no_full_auc_replicates",
              "undefined_full_auc_times", "undefined_observed_area_times", "incomplete_path_times",
              "complete_path_undefined_auc_times")
    with (root / "roc-missingness.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("setting", "method", "target", *fields))
        for (setting, method, target), counts in sorted(roc_totals.items()):
            assert counts["completed_datasets"] + counts["incomplete_datasets"] == counts["intended_replicates"]
            assert (counts["complete_full_auc_replicates"] + counts["partial_full_auc_replicates"] +
                    counts["no_full_auc_replicates"]) == counts["completed_datasets"]
            writer.writerow((setting, method, target, *(counts[f] for f in fields)))
    print("Selection and ROC missingness saved separately from candidate solver failures")


if __name__ == "__main__":
    main(sys.argv[1])
