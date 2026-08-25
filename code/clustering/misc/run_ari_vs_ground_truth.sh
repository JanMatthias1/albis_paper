#!/bin/bash
#SBATCH --job-name=sim_ari_vs_ground_truth
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/sim_ari_vs_ground_truth_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=08:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared      # uncomment / edit to match your cluster's partition name

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs

source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

MODALITIES=(spot bin cell)
MODALITY="${MODALITIES[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "Modality: ${MODALITY}"
echo "Command: ${PYTHON_BIN} sim_paper/code/clustering/ari_vs_ground_truth.py --modality ${MODALITY} $*"
"${PYTHON_BIN}" sim_paper/code/clustering/ari_vs_ground_truth.py --modality "${MODALITY}" "$@"

echo "Job finished: $(date)"
