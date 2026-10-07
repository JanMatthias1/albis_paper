#!/bin/bash
#SBATCH --job-name=real_data_qc
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/real_data_qc/logs/real_data_qc_%A_%a.out
#SBATCH --array=0-1
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/albis_paper/code/real_data_qc/logs

source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/real_data_qc/_env.sh

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

SLICES=(non_diseased_lung lung_cancer)
SLICE="${SLICES[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "Slice: ${SLICE}"
echo "Command: ${PYTHON_BIN} albis_paper/code/real_data_qc/xenium_qc.py --slice ${SLICE} $*"
"${PYTHON_BIN}" albis_paper/code/real_data_qc/xenium_qc.py --slice "${SLICE}" "$@"

echo "Job finished: $(date)"
