#!/bin/bash
#SBATCH --job-name=sim_count_distribution
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/sim_count_distribution_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs

source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

MODALITIES=(spot bin cell)
MODALITY="${MODALITIES[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "Modality: ${MODALITY}"
echo "Command: ${PYTHON_BIN} sim_paper/code/count_distribution/count_distribution.py --modality ${MODALITY} $*"
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py --modality "${MODALITY}" "$@"

echo "Job finished: $(date)"
