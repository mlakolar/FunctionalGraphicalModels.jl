import csv
import pathlib
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root=pathlib.Path(sys.argv[1])
with (root/"paired.csv").open(newline="") as f: rows=list(csv.DictReader(f))
x=[float(r["default_frobenius"]) for r in rows]
y=[float(r["selected_frobenius"]) for r in rows]
fig,ax=plt.subplots(figsize=(5.5,5),layout="constrained")
lo,hi=min(x+y)*.75,max(x+y)*1.3
ax.plot([lo,hi],[lo,hi],color="gray",linestyle="--",linewidth=1)
ax.scatter(x,y,s=30,alpha=.8)
ax.set(xscale="log",yscale="log",xlim=(lo,hi),ylim=(lo,hi),
    xlabel="Default scaled lasso: normalized squared Frobenius loss",
    ylabel="Held-out tuned scaled lasso: same loss",
    title="Scaled-lasso calibration · 48 paired datasets\nSix candidate penalties; no latent truth used for selection")
ax.grid(alpha=.2)
for ext in ("png","svg"): fig.savefig(root/f"calibration.{ext}",dpi=180)
