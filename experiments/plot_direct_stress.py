import csv
import pathlib
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = pathlib.Path(sys.argv[1])
with (root / "summary.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))
fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), sharey=True, layout="constrained")
for ax, (curvature, noise) in zip(axes.flat, ((0, .5), (0, 2), (4, .5), (4, 2))):
    for method, label in (("aclime", "ACLIME"), ("l2_aclime", "ℓ₂-ACLIME"), ("scaled_lasso", "Scaled lasso")):
        chosen = [r for r in rows if r["method"] == method and r["quantity"] == "normalized_squared_frobenius"
            and float(r["curvature"]) == curvature and float(r["noise_variance"]) == noise]
        chosen.sort(key=lambda r: (int(r["degree"]), int(r["visits"])))
        ax.errorbar(range(len(chosen)), [float(r["mean"]) for r in chosen],
            yerr=[float(r["mcse"]) for r in chosen], marker="o", capsize=3, label=label)
    ax.set_xticks(range(4), ("d=2, m=8", "d=2, m=30", "d=10, m=8", "d=10, m=30"), rotation=15)
    ax.set_title(f"Curvature {curvature}; noise variance {noise:g}", fontsize=11)
    ax.tick_params(labelsize=9)
    ax.set_yscale("log")
    ax.grid(alpha=.2)
handles, labels = axes[0,0].get_legend_handles_labels()
fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False)
fig.supylabel("Normalized squared Frobenius loss (log scale)", fontsize=11)
fig.suptitle("Known-target controls · n=200, p=50\nThree paired replicates; bars show estimated Monte Carlo standard errors", fontsize=12)
for extension in ("png", "svg"):
    fig.savefig(root / f"selected-losses.{extension}", dpi=180)
