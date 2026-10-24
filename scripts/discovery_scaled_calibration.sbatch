#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail
module load parallel/20240522
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export DFGM_SOURCE_DIR=${1:?source required}
export DFGM_OUTPUT_DIR=${2:?output required}
mkdir -p "$DFGM_OUTPUT_DIR"
# Validate the full new path once before dispatching the remaining cases.
srun --exclusive --exact -N1 -n1 -c1 --mem=1800M \
  bash "$DFGM_SOURCE_DIR/scripts/scaled_calibration_worker.sh" 1
worker_status=0
parallel --jobs 2 --halt never --joblog "$DFGM_OUTPUT_DIR/parallel.tsv" \
 'srun --exclusive --exact -N1 -n1 -c1 --mem=1800M bash "$DFGM_SOURCE_DIR/scripts/scaled_calibration_worker.sh" {}' \
 ::: {2..16} || worker_status=$?
python3 "$DFGM_SOURCE_DIR/experiments/summarize_scaled_calibration.py" "$DFGM_OUTPUT_DIR" "${3:?original stress outputs required}"
export MPLCONFIGDIR="$DFGM_OUTPUT_DIR/.mpl-cache"
/home1/mkolar/juliaPkg/runs/20261024-fgm-development/plot-env/bin/python \
  "$DFGM_SOURCE_DIR/experiments/plot_scaled_calibration.py" "$DFGM_OUTPUT_DIR"
exit "$worker_status"
