#!/bin/bash
# Orchestration only: every scientific worker is an exclusive Slurm job step.
set -euo pipefail
module load parallel/20240522
job_id=${1:?existing owned allocation required}
export DFGM_SOURCE_DIR=${2:?source required}
export DFGM_OUTPUT_DIR=${3:?output required}
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$DFGM_OUTPUT_DIR"
export DFGM_ALLOCATION_ID="$job_id"
parallel --jobs 2 --halt never --joblog "$DFGM_OUTPUT_DIR/parallel.tsv" \
  'srun --jobid="$DFGM_ALLOCATION_ID" --exclusive --exact -N1 -n1 -c1 --mem=1800M --time=00:03:00 bash "$DFGM_SOURCE_DIR/scripts/direct_stress_worker.sh" {}' \
  ::: {1..16}
