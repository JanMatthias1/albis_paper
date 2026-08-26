#!/bin/bash
#SBATCH --job-name=sim_count_distribution_compare
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/sim_count_distribution_compare_%A_%a.out
#SBATCH --array=0-1
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

SLICES=(non_diseased_lung lung_cancer)
SLICE="${SLICES[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "Compare slice: ${SLICE}"
echo "Command: ${PYTHON_BIN} sim_paper/code/count_distribution/count_distribution.py --modality cell --compare-input sim_paper/data/real_data_qc/${SLICE}/${SLICE}_qc.h5ad --compare-label ${SLICE} $*"
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality cell \
    --compare-input "sim_paper/data/real_data_qc/${SLICE}/${SLICE}_qc.h5ad" \
    --compare-label "${SLICE}" \
    "$@"

echo "Job finished: $(date)"
