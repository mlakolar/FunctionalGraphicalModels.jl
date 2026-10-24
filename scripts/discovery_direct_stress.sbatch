#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
set -euo pipefail
module load parallel/20240522
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
source_dir=${1:?source directory required}
output_dir=${2:?output directory required}
mkdir -p "$output_dir"
cd "$source_dir"
export DFGM_SOURCE_DIR="$source_dir" DFGM_OUTPUT_DIR="$output_dir"
worker_status=0
parallel --jobs "$SLURM_NTASKS" --joblog "$output_dir/parallel.tsv" --halt never \
  'srun --exclusive --exact -N1 -n1 -c1 --mem=1800M bash scripts/direct_stress_worker.sh {}' \
  ::: {1..16} || worker_status=$?
if [[ -n ${DFGM_PILOT_OUTPUT:-} ]]; then
    /home1/mkolar/juliaPkg/runs/20261024-fgm-development/plot-env/bin/python \
      experiments/plot_dfgm.py "$DFGM_PILOT_OUTPUT"
fi
exit "$worker_status"
