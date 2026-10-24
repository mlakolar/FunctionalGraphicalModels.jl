#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:02:00
set -euo pipefail
root=${1:?run root required}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR="$root/outputs/.mpl-cache"
python=/home1/mkolar/juliaPkg/runs/20261024-fgm-development/plot-env/bin/python
"$python" "$root/analysis-source/plot_direct_stress_v2.py" "$root/outputs/stress"
"$python" "$root/analysis-source/plot_dfgm_v2.py" "$root/outputs/pilot"
