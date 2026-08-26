#!/bin/bash
#SBATCH --job-name=sim_count_distribution_compare_noisy_visium
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/sim_count_distribution_compare_noisy_visium_%j.out
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

# Optional: OUT_TAG (e.g. exported via `sbatch --export=ALL,OUT_TAG=packing_pf0p06 ...`) reads
# the simulation from data/noisy/<tag>/ and writes comparison plots into a <tag> subfolder
# instead of overwriting the top-level simulation_spot_z.h5ad / spot_noisy_vs_breast_cancer_visium/*.png.
OUT_TAG="${OUT_TAG:-}"
if [[ -n "${OUT_TAG}" ]]; then
    INPUT_PATH="sim_paper/data/noisy/${OUT_TAG}/simulation_spot_z.h5ad"
    OUTPUT_DIR="sim_paper/data/count_distribution/spot_noisy_vs_breast_cancer_visium/${OUT_TAG}"
else
    INPUT_PATH="sim_paper/data/noisy/simulation_spot_z.h5ad"
    OUTPUT_DIR="sim_paper/data/count_distribution/spot_noisy_vs_breast_cancer_visium"
fi

echo "Command: ${PYTHON_BIN} sim_paper/code/count_distribution/count_distribution.py --modality spot --input ${INPUT_PATH} --compare-input sim_paper/data/real_data_qc/breast_cancer_visium/breast_cancer_visium_qc.h5ad --compare-label breast_cancer_visium --output-dir ${OUTPUT_DIR} $*"
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality spot \
    --input "${INPUT_PATH}" \
    --compare-input "sim_paper/data/real_data_qc/breast_cancer_visium/breast_cancer_visium_qc.h5ad" \
    --compare-label "breast_cancer_visium" \
    --output-dir "${OUTPUT_DIR}" \
    "$@"

echo "Job finished: $(date)"
