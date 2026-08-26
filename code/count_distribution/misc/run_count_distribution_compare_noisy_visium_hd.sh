#!/bin/bash
#SBATCH --job-name=sim_count_distribution_compare_noisy_visium_hd
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/sim_count_distribution_compare_noisy_visium_hd_%A_%a.out
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

SAMPLES=(breast_cancer_visium_hd human_pancreas_visium_hd)
SAMPLE="${SAMPLES[${SLURM_ARRAY_TASK_ID:-0}]}"

# Optional: OUT_TAG (e.g. exported via `sbatch --export=ALL,OUT_TAG=theta_200 ...`) reads
# the simulation from data/noisy/<tag>/ and writes comparison plots into a <tag> subfolder
# instead of overwriting the top-level simulation_bin_z.h5ad / bin_noisy_vs_<sample>/*.png.
OUT_TAG="${OUT_TAG:-}"
if [[ -n "${OUT_TAG}" ]]; then
    INPUT_PATH="sim_paper/data/noisy/${OUT_TAG}/simulation_bin_z.h5ad"
    OUTPUT_DIR="sim_paper/data/count_distribution/bin_noisy_vs_${SAMPLE}/${OUT_TAG}"
else
    INPUT_PATH="sim_paper/data/noisy/simulation_bin_z.h5ad"
    OUTPUT_DIR="sim_paper/data/count_distribution/bin_noisy_vs_${SAMPLE}"
fi

echo "Sample: ${SAMPLE}"
echo "Command: ${PYTHON_BIN} sim_paper/code/count_distribution/count_distribution.py --modality bin --input ${INPUT_PATH} --compare-input sim_paper/data/real_data_qc/${SAMPLE}/${SAMPLE}_qc.h5ad --compare-label ${SAMPLE} --output-dir ${OUTPUT_DIR} $*"
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality bin \
    --input "${INPUT_PATH}" \
    --compare-input "sim_paper/data/real_data_qc/${SAMPLE}/${SAMPLE}_qc.h5ad" \
    --compare-label "${SAMPLE}" \
    --output-dir "${OUTPUT_DIR}" \
    "$@"

echo "Job finished: $(date)"
