#!/bin/bash
#SBATCH --job-name=crossmod_gene_metrics
#SBATCH --partition=shared
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/cross_modality_alignment/independent_offsets_shift3x/gene_metrics_%j.out
# Run directly with bash, or submit with sbatch. Extra arguments go to Python.
set -euo pipefail
ROOT=/dcs04/hicks/data/Jan/sim_project
CODE="$ROOT/sim_paper/code/applications_albis/cross_modality_alignment"
export MPLCONFIGDIR="/tmp/crossmod-gene-mpl-${SLURM_JOB_ID:-local}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
"$ROOT/sim_paper/env/albis-tutorial/bin/python" -u "$CODE/evaluation_gene_metrics.py" "$@"
