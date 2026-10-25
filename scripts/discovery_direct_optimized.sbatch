#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${1:?optimized source required}"
cp "${DFGM_BACKEND_SOURCE:?resolved backend source required}/Manifest.toml" Manifest.toml
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
"$JULIA" --project -e 'using FunctionalGraphicalModels, LinearAlgebra, Random, Test; import FunctionalDataAnalysis as FDA; include("test/direct.jl")'
"$JULIA" --project experiments/benchmark_dfgm_common.jl --out "${2:?output directory required}"
