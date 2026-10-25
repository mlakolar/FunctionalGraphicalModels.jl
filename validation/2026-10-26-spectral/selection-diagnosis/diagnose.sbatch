#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
root=/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct
python3 "$root/spectral-selection-diagnosis-v1/diagnose.py" "$root" "$root/outputs/spectral-selection-diagnosis-$SLURM_JOB_ID"
