#!/bin/bash
#SBATCH --job-name=brain_xenium_qc
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/real_data_qc/logs/brain_xenium_qc_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
# QC for the three Xenium human brain slices (10x v1.3.0, FFPE, with add-on panel), same
# 4-MAD QC as the Figure 2 Xenium lung references. Raw: data/real_data/brain_samples_xenium/
# Output: data/real_data_qc/brain_samples_xenium/<slice>/

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/albis_paper/code/real_data_qc/logs

source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/real_data_qc/_env.sh

cd /dcs04/hicks/data/Jan/sim_project

SLICES=(brain_healthy brain_glioblastoma brain_alzheimers)
SLICE="${SLICES[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "Job started: $(date)  Host: $(hostname)  Slice: ${SLICE}"
"${PYTHON_BIN}" albis_paper/code/real_data_qc/xenium_qc.py --slice "${SLICE}" \
  --output-dir "albis_paper/data/real_data_qc/brain_samples_xenium/${SLICE}" "$@"
echo "Job finished: $(date)"
