#!/bin/bash
#SBATCH --job-name=sim_count_distribution_compare_noisy
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/sim_count_distribution_compare_noisy_%A_%a.out
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

# Optional: OUT_TAG (e.g. exported via `sbatch --export=ALL,OUT_TAG=log_mu_-2.5 ...`) reads
# the simulation from data/noisy/<tag>/ and writes comparison plots into a <tag> subfolder
# instead of overwriting the top-level simulation_cell_z.h5ad / cell_noisy_vs_<slice>/*.png.
OUT_TAG="${OUT_TAG:-}"
if [[ -n "${OUT_TAG}" ]]; then
    INPUT_PATH="sim_paper/data/noisy/${OUT_TAG}/simulation_cell_z.h5ad"
    OUTPUT_DIR="sim_paper/data/count_distribution/cell_noisy_vs_${SLICE}/${OUT_TAG}"
else
    INPUT_PATH="sim_paper/data/noisy/simulation_cell_z.h5ad"
    OUTPUT_DIR="sim_paper/data/count_distribution/cell_noisy_vs_${SLICE}"
fi

echo "Compare slice: ${SLICE}"
echo "Command: ${PYTHON_BIN} sim_paper/code/count_distribution/count_distribution.py --modality cell --input ${INPUT_PATH} --compare-input sim_paper/data/real_data_qc/${SLICE}/${SLICE}_qc.h5ad --compare-label ${SLICE} --output-dir ${OUTPUT_DIR} $*"
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality cell \
    --input "${INPUT_PATH}" \
    --compare-input "sim_paper/data/real_data_qc/${SLICE}/${SLICE}_qc.h5ad" \
    --compare-label "${SLICE}" \
    --output-dir "${OUTPUT_DIR}" \
    "$@"

echo "Job finished: $(date)"
