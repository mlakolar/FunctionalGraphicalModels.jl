"""Static paper-draft pilot panels from validated reporting CSVs (matplotlib).
Run under Slurm, after summarize_dfgm.py. Gaps remain missing; no imputation.
"""
import csv
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = pathlib.Path(sys.argv[1])
paper = "--paper" in sys.argv[2:]
with (root / "summary.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))
methods = ("aclime", "l2_aclime", "scaled_lasso", "DFGM_C", "DFGM_G")
labels = ("Direct ACLIME", "Direct ℓ₂-ACLIME", "Direct scaled lasso", "DFGM-C", "DFGM-G")
if any(r["method"] == "scaled_lasso_tuned" for r in rows):
    methods += ("scaled_lasso_tuned",)
    labels += ("Scaled lasso, BIC menu",)
if paper:
    methods += ("DGM_C", "DGM_G", "GGL", "Naive_C", "Naive_G")
    labels += ("DGM-C", "DGM-G", "GGL", "Naive-C", "Naive-G")
for target in ("nominal", "actual"):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8) if paper else (9, 6), sharey=True, layout="constrained")
    for ax, indices, title in zip(axes.flat, ((1,2,3),(4,5,6),(7,8,9),(10,11,12)),
            ("Common, p=50", "Common, p=100", "Irregular, p=50", "Irregular, p=100")):
        for method, label in zip(methods, labels):
            selected = {int(r["setting"]): r for r in rows
                if r["method"] == method and r["target"] == target and r["quantity"] == "operator"}
            points = [float(selected[i]["mean"]) if i in selected else float("nan") for i in indices]
            errors = [float(selected[i]["mcse"]) if i in selected else float("nan") for i in indices]
            ax.errorbar(range(3), points, yerr=errors, marker="o", linewidth=1.4, label=label)
        ax.set_title(title)
        ax.set_xticks(range(3), ("10", "25", "50") if indices[0]<7 else ("6–9", "20–30", "50"))
        ax.set_xlabel("Visits per curve")
        ax.set_ylabel("Mean operator norm loss")
        ax.grid(alpha=.2)
    handles, names = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, names, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle(f"DFGM comparison · {target} target\nMeans and Monte Carlo SE; available replicate counts in summary.csv" if paper else
                f"Exploratory DFGM pilot · {target} target\nOne realization; initial calibration; gaps mark unavailable complete-time summaries")
    for extension in ("png", "svg"):
        fig.savefig(root / f"operator-{target}.{extension}", dpi=180)
    plt.close(fig)
