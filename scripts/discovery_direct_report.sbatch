#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:03:00
set -euo pipefail
root=${1:?run root required}
source_dir=${2:?source required}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MPLCONFIGDIR="$root/outputs/.mpl-cache"
python3 "$source_dir/experiments/summarize_dfgm.py" "$root/outputs/pilot"
python3 "$root/analysis-source/summarize_direct_stress.py" "$root/outputs/stress"
/home1/mkolar/juliaPkg/runs/20261024-fgm-development/plot-env/bin/python \
  "$source_dir/experiments/plot_dfgm.py" "$root/outputs/pilot"
# Verify the reporting inputs and expected candidate/time coverage without refits.
python3 - "$root" <<'PY'
import csv, json, math, pathlib, sys, tomllib
r=pathlib.Path(sys.argv[1]); checked=[]
for index in range(1,13):
    d=r/"outputs"/"pilot"/f"setting-{index}-rep-1"
    if not (d/"complete.toml").exists():
        checked.append(dict(setting=index,status="incomplete")); continue
    settings=tomllib.loads((d/"settings.toml").read_text())
    receipt=tomllib.loads((d/"complete.toml").read_text())
    with (d/"fits.csv").open(newline="") as f: fits=list(csv.DictReader(f))
    expected=28 if settings["design"]=="common" else 13
    assert len(fits)==expected,(index,len(fits),expected)
    assert sum(x["status"]=="failed" for x in fits)==receipt["fit_failures"]
    with (d/"metrics.csv").open(newline="") as f: metrics=list(csv.DictReader(f))
    q=len(settings["evaluation_grid"])
    assert len(metrics)==2*q*sum(x["status"]=="ok" for x in fits)
    assert all(math.isfinite(float(x[k])) for x in metrics for k in
               ("operator","frobenius","normalized_squared_frobenius"))
    checked.append(dict(setting=index,status=receipt["status"],candidate_fits=len(fits),
                        time_points=q,metric_rows=len(metrics),fit_failures=receipt["fit_failures"]))
(r/"outputs"/"report-validation.json").write_text(json.dumps(dict(
    pilot=checked,checks="expected candidate/time counts; recorded failures; finite precision losses; stress coverage and selected-fit uniqueness; no numerical refits"),indent=2)+"\n")
print(json.dumps(checked))
PY
