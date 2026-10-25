#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513|epyc-7542
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
root=/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct
python3 "$root/spectral-report-full-v3/report.py" "$root" "$root/outputs/spectral-full-current-v3-$SLURM_JOB_ID"
