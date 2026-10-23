#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
# Submit from the package directory after checking these project-specific settings.
set -euo pipefail
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=2
export JULIA_NUM_PRECOMPILE_TASKS=2
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
FGM_JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
cd "$SLURM_SUBMIT_DIR"
"$FGM_JULIA" scripts/instantiate.jl
"$FGM_JULIA" --project -e 'using Pkg; Pkg.test()'
"$FGM_JULIA" --project experiments/qiao2019.jl --smoke --model all --out "results/qiao-smoke-$SLURM_JOB_ID"
"$FGM_JULIA" --project experiments/other_methods_smoke.jl "results/other-methods-$SLURM_JOB_ID"
