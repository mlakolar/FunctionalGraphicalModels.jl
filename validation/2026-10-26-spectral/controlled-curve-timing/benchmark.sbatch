#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513|epyc-7542
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
set -euo pipefail
root=/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=4 JULIA_NUM_GC_THREADS=1 JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_HEAP_SIZE_HINT=5000M
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia --project="$root/common-fast-v2-source" "$root/spectral-curve-timing-v1/benchmark.jl" "$root/outputs/spectral-curve-timing-$SLURM_JOB_ID"
