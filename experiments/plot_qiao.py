#!/usr/bin/env python3
"""Plot observed Qiao comparison paths; requires Python 3.11 and matplotlib."""
import argparse
import csv
from pathlib import Path
import tomllib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics", type=Path)
    parser.add_argument("output", type=Path, help="new output stem (writes .svg and .png)")
    parser.add_argument("--model", default="ar2")
    parser.add_argument("--replicate", type=int, default=1)
    args = parser.parse_args()
    outputs = [args.output.with_suffix(ext) for ext in (".svg", ".png")]
    if any(path.exists() for path in outputs):
        parser.error("output already exists")
    with args.metrics.open(newline="") as stream:
        rows = [row for row in csv.DictReader(stream)
                if row["model"] == args.model and int(row["replicate"]) == args.replicate]
    if not rows:
        parser.error("no rows for the requested realization")
    with args.metrics.with_name("settings.toml").open("rb") as stream:
        settings = tomllib.load(stream)
    fig, ax = plt.subplots(figsize=(8, 6.2), constrained_layout=True)
    colors = ["#111111", "#D55E00", "#009E73", "#0072B2", "#CC79A7"]
    styles = ["-", "--", ":", "-.", "--"]
    for method, color, style in zip(("fglasso", "ALL", "NEVER", "HALF", "PCA"), colors, styles):
        attempted = [row for row in rows if row["method"] == method]
        valid = [row for row in attempted if row["valid"].lower() == "true"]
        valid.sort(key=lambda row: (float(row["fpr"]), int(row["lambda_index"])))
        label = f"{method} ({len(valid)}/{len(attempted)} valid)"
        ax.plot([float(row["fpr"]) for row in valid], [float(row["tpr"]) for row in valid],
                marker="o", markersize=4, linewidth=1.7, color=color, linestyle=style, label=label)
    ax.plot([0, 1], [0, 1], color="0.75", linestyle=":", linewidth=1, label="Chance reference")
    ax.set(xlim=(-0.015, 1.015), ylim=(-0.015, 1.015),
           xlabel="False positive rate", ylabel="True positive rate")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.18)
    ax.legend(loc="lower right", fontsize=9)
    fig.suptitle(f"Qiao {args.model} pilot: n={settings['n']}, p={settings['p']}, T={settings['ntime']}", fontsize=13)
    ax.set_title(f"One dataset; seed {rows[0]['seed']}; observed fitted points only", fontsize=10, pad=10)
    fig.text(0.5, 0.005, "Coarse regularization paths. No synthetic endpoints or published-table replication.",
             ha="center", fontsize=8, color="0.35")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for path in outputs:
        fig.savefig(path, dpi=180, bbox_inches="tight")
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
