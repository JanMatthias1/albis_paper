#!/bin/bash
#SBATCH --job-name=fig2_bin_vs_breast_cancer_visium_hd_with_batch
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/fig2_bin_vs_breast_cancer_visium_hd_with_batch_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Same pairing as bin_vs_breast_cancer_visium_hd.sh (bin vs. Visium HD
# breast_cancer_visium_hd), but with --use-batch-effect: real data has no
# "pre-batch" version, whatever technical noise it carries is just baked
# into its counts, so stripping the sim's synthetic batch effect (the
# default) makes sim artificially cleaner than real can ever be. Now that
# --slice-id restricts to a single slice (5, see below), this is the fair
# apples-to-apples comparison. Writes to figure_2_with_batch/ instead of
# figure_2/ so both variants are kept side by side. See
# bin_vs_breast_cancer_visium_hd.sh for the sim config rationale.
#
# --slice-id 5: restricts the sim (10 z-plane slices through the sphere) to
# one middle slice, since real data is a single tissue section. Slices near
# the sphere's poles (0 and 9) sample a much smaller cross-section and are
# degenerate (up to 20x higher empty-bin fraction) -- pooling all 10 masked
# this. Slice 5 (near-equator) is representative.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_log_mu_0.0"
MODALITY="bin"
REAL_LABEL="breast_cancer_visium_hd"
SIM_RAW="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/figure_2_with_batch/${MODALITY}_vs_${REAL_LABEL}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um 2050 \
        --base-gene-lognormal 0.0 0.7 \
        --out-tag "${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
fi

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" --slice-id 5 --use-batch-effect \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/full_panel"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" --slice-id 5 --use-batch-effect \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/hvg_matched"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" --slice-id 5 --use-batch-effect \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/qc_filtered"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" --slice-id 5 --use-batch-effect \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

echo "[done] All 4 modes written under ${OUT_ROOT}"
