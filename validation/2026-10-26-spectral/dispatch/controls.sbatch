#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513|epyc-7542
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
root=/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct
module load parallel/20240522
parallel --line-buffer --jobs 4 --halt never --joblog "$root/outputs/spectral-v5/dispatch/control-joblog-$SLURM_JOB_ID.tsv" \
 srun --exclusive --exact -N1 -n1 -c1 --mem=3500M bash "$root/spectral-dispatch-v1/control-worker.sh" '{1}' :::: "$1"
