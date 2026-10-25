"""Validate complete coverage, observed-data selection and fitted ROC accounting.

Run under Slurm. Optional source locks enforce the production protocol. Optional
independent CSVs come from saved precision/covariance arrays, without refits.
"""
import argparse
import csv
import json
import math
import pathlib
import sys
import tomllib
from collections import Counter, defaultdict


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def near(a, b, *, atol=1e-8, rtol=1e-9):
    a, b = float(a), float(b)
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if math.isinf(a) or math.isinf(b):
        return a == b
    return math.isclose(a, b, abs_tol=atol, rel_tol=rtol)


def source_settings(settings, locks):
    variants = [v for v in locks["variants"]
                if settings["source_sha256"] == v["source_sha256"]]
    assert len(variants) == 1, "unrecognized or mixed numerical source hashes"
    variant = variants[0]
    assert settings["dependency_manifest_sha256"] == locks["dependency_manifest_sha256"]
    assert settings["julia_version"] == locks["julia_version"]
    setting, replicate = settings["setting_index"], settings["replicate"]
    assert 1 <= setting <= 12 and 1 <= replicate <= 100
    assert setting in variant["settings"], "source variant assigned to unexpected setting"
    assert settings["seed"] == 2026102500 + 10000 * setting + replicate
    assert settings["n"] == 200
    assert settings["p"] == (50 if ((setting - 1) // 3) % 2 == 0 else 100)
    assert settings["design"] == ("common" if setting <= 6 else "irregular")
    schedules = ("T10", "T25", "T50") if setting <= 6 else ("T6_9", "T20_30", "T50")
    assert settings["schedule"] == schedules[(setting - 1) % 3]
    assert settings["mode"] == "paper" and settings["integration_points"] == 101
    assert settings["work_points"] == 31 and settings["noise_variance"] == .5
    assert settings["path_reuse"] and settings["threads_per_worker"] == 4
    expected_fast = setting <= 6 if variant.get("common_fast_by_design") else variant["common_fast"]
    assert settings.get("common_fast", False) == expected_fast
    assert settings["bias_constant_menu"] == [.01, .1, 1., 10., 20.]
    assert len(settings["K2_menu"]) == 31
    assert near(settings["K2_menu"][0], .03) and near(settings["K2_menu"][-1], 2.)
    assert settings["bandwidths"] == [.15, .25, .4]
    assert settings["threshold"] == 1e-6
    assert settings["K1"] == 2. and settings["lp_solver"] == "highs"
    q = len(settings["evaluation_grid"])
    assert q == ((10, 25, 50)[(setting - 1) % 3] if setting <= 6 else 21)
    counts = {"aclime": 155, "l2_aclime": 155, "scaled_lasso": 1,
              "scaled_lasso_tuned": 6, "DFGM_C": 31, "DFGM_G": 31}
    if setting <= 6:
        counts.update({m: 31 for m in ("DGM_C", "DGM_G", "GGL", "Naive_C", "Naive_G")})
    subset = settings.get("method_subset", "all")
    assert subset == variant.get("method_subset", "all")
    if subset == "scaled_original":
        assert replicate > 20, "retained first20 must contain both proposals"
        counts.pop("aclime")
        counts.pop("l2_aclime")
    assert settings["expected_candidates"] == counts
    return variant["name"]


def roc_area(points):
    usable = [(x, y) for x, y in points if math.isfinite(x) and math.isfinite(y)]
    if not usable:
        return (math.nan,) * 4
    # Multiple fitted points at one FPR span a vertical segment of zero area.
    highest = {}
    for x, y in usable:
        highest[x] = max(highest.get(x, -math.inf), y)
    ordered = sorted(highest.items())
    area = sum((x1-x0) * (y1+y0)/2 for (x0, y0), (x1, y1) in zip(ordered, ordered[1:]))
    full = area if (0., 0.) in usable and (1., 1.) in usable else math.nan
    return area, full, ordered[0][0], ordered[-1][0]


def validate(root, *, allow_recorded_fit_failures=False, source_locks=None, independent_root=None):
    root = pathlib.Path(root)
    locks = json.loads(pathlib.Path(source_locks).read_text()) if source_locks else None
    tasks = [tuple(map(int, line.split("\t")[:2])) for line in (root / "tasks.tsv").read_text().splitlines() if line]
    assert len(tasks) == len(set(tasks)), "duplicate intended dataset"
    if locks:
        assert set(tasks) == {(s, r) for s in range(1, 13) for r in range(1, 101)}
    checks = []
    for setting, replicate in tasks:
        folder = root / f"setting-{setting}-rep-{replicate}"
        if not (folder / "complete.toml").exists():
            checks.append(dict(setting=setting, replicate=replicate, status="incomplete"))
            continue
        settings = tomllib.loads((folder / "settings.toml").read_text())
        receipt = tomllib.loads((folder / "complete.toml").read_text())
        assert settings["setting_index"] == setting and settings["replicate"] == replicate
        variant = source_settings(settings, locks) if locks else "not_checked"
        counts = settings["expected_candidates"]
        fits = read_csv(folder / "fits.csv")
        keyed = {(r["method"], int(r["candidate"])): r for r in fits}
        expected = {(m, c) for m, n in counts.items() for c in range(1, n+1)}
        assert len(fits) == len(keyed) and set(keyed) == expected, (folder, "fit coverage")
        assert all(r["status"] in ("ok", "failed") for r in fits)
        failures = sum(r["status"] == "failed" for r in fits)
        assert failures == receipt["fit_failures"]
        assert receipt["status"] == ("completed_with_failures" if failures else "completed")
        successful = {key: r for key, r in keyed.items() if r["status"] == "ok"}
        q, p = len(settings["evaluation_grid"]), settings["p"]
        metrics = read_csv(folder / "metrics.csv")
        indexed = {(r["method"], r["target"], int(r["time_index"]), int(r["candidate"])): r for r in metrics}
        expected_metrics = {(m, target, t, c) for m, c in successful
                            for target in ("actual", "nominal") for t in range(1, q+1)}
        assert len(metrics) == len(indexed) and set(indexed) == expected_metrics, (folder, "metric coverage")
        grouped = defaultdict(list)
        for row in metrics:
            assert all(math.isfinite(float(row[k])) for k in ("operator", "frobenius", "normalized_squared_frobenius"))
            assert int(row["edges"]) + int(row["nonedges"]) == p*(p-1)//2
            for quantity, denominator in (("tpr", "edges"), ("fpr", "nonedges")):
                value = float(row[quantity])
                assert (math.isnan(value) if int(row[denominator]) == 0 else math.isfinite(value) and 0 <= value <= 1)
            grouped[row["method"], row["target"], int(row["time_index"])].append(row)
        with (folder / "selected.csv").open(newline="") as stream:
            selected_rows = list(csv.reader(stream))
        selections = {(m, int(t)): (int(c), status) for m, t, c, status in selected_rows}
        successful_methods = {method for method, _ in successful}
        expected_selections = {(m, t) for m in successful_methods for t in range(1, q+1)}
        assert len(selected_rows) == len(selections) and set(selections) == expected_selections
        missing = Counter({m: q for m in counts})
        for (method, t), (candidate, status) in selections.items():
            points = grouped[method, "actual", t]
            best = min(points, key=lambda r: float(r["bic"]))  # first minimum, original path order
            if math.isfinite(float(best["bic"])):
                assert status == "selected" and candidate == int(best["candidate"]), (folder, method, t, "BIC argmin")
                missing[method] -= 1
            else:
                assert candidate == 0 and status == "no_positive_definite_candidate"
            for row in points:
                other = indexed[method, "nominal", t, int(row["candidate"])]
                assert float(row["bic"]) == float(other["bic"]), "target-dependent selection criterion"
        auc = read_csv(folder / "auc.csv")
        auc_keys = {(r["method"], r["target"], int(r["time_index"])) for r in auc}
        assert len(auc) == len(auc_keys) and auc_keys == {(m, target, t) for m in counts for target in ("actual", "nominal") for t in range(1, q+1)}
        for row in auc:
            method, target, t = row["method"], row["target"], int(row["time_index"])
            path = grouped[method, target, t]
            n_expected = counts[method]
            if method in ("aclime", "l2_aclime") and len(settings["bias_constant_menu"]) > 1:
                candidate, status = selections.get((method, t), (0, "missing"))
                chosen = float(successful[method, candidate]["bias_constant"]) if status == "selected" else math.nan
                assert near(row["selected_bias_constant"], chosen)
                path = [r for r in path if float(successful[method, int(r["candidate"])]["bias_constant"]) == chosen]
                n_expected = len(settings["K2_menu"])
            assert int(row["expected_candidates"]) == n_expected
            assert int(row["successful_candidates"]) == len(path)
            area, full, lo, hi = roc_area([(float(r["fpr"]), float(r["tpr"])) for r in path])
            if len(path) != n_expected:
                full = math.nan
            assert all(near(row[k], value) for k, value in
                       (("observed_area", area), ("full_auc", full), ("minimum_fpr", lo), ("maximum_fpr", hi))), (folder, method, t, "ROC accounting")
        assert (folder / "covariances.bin").stat().st_size > 0
        restarts = read_csv(folder / "numerical-restarts.csv") if (folder / "numerical-restarts.csv").exists() else []
        if "successful_cold_column_restarts" in receipt:
            assert sum(int(r["column_restarts"]) for r in restarts) == receipt["successful_cold_column_restarts"]
            assert all("ALMOST_OPTIMAL" in r["initial_error"] and r["recovery"] == "independent_cold_OPTIMAL_same_tolerances" for r in restarts)
        independent_count = 0
        max_bic_error = max_loss_error = 0.
        if independent_root:
            independent = pathlib.Path(independent_root) / folder.name / "independent-metrics.csv"
            if independent.exists():
                recomputed = read_csv(independent)
                keys = {(r["method"], r["target"], int(r["time_index"]), int(r["candidate"])) for r in recomputed}
                assert len(recomputed) == len(keys) and keys == expected_metrics
                for row in recomputed:
                    actual = indexed[row["method"], row["target"], int(row["time_index"]), int(row["candidate"])]
                    for key in ("bic", "operator", "frobenius", "normalized_squared_frobenius", "fpr", "tpr", "edges", "nonedges"):
                        assert near(row[key], actual[key]), (folder, key, row)
                        if math.isfinite(float(row[key])) and key in ("bic", "operator", "frobenius", "normalized_squared_frobenius"):
                            error = abs(float(row[key]) - float(actual[key]))
                            if key == "bic": max_bic_error = max(max_bic_error, error)
                            else: max_loss_error = max(max_loss_error, error)
                independent_count = len(recomputed)
        checks.append(dict(setting=setting, replicate=replicate, status=receipt["status"], source_variant=variant,
                           candidate_fits=len(fits), fit_failures=failures, metric_rows=len(metrics),
                           missing_selected_times=dict(missing), independent_metric_rows=independent_count,
                           max_independent_bic_error=max_bic_error, max_independent_loss_error=max_loss_error,
                           successful_cold_column_restarts=receipt.get("successful_cold_column_restarts", 0)))
    (root / "validation.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks))
    return all(r["status"] in ("completed", "completed_with_failures") and
               (allow_recorded_fit_failures or r["fit_failures"] == 0) for r in checks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root")
    parser.add_argument("--allow-recorded-fit-failures", action="store_true")
    parser.add_argument("--source-locks")
    parser.add_argument("--independent-root")
    args = parser.parse_args()
    sys.exit(0 if validate(args.root, allow_recorded_fit_failures=args.allow_recorded_fit_failures,
                          source_locks=args.source_locks, independent_root=args.independent_root) else 1)
