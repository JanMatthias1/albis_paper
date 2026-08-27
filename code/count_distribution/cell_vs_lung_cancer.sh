#!/bin/bash
#SBATCH --job-name=fig2_cell_vs_lung_cancer
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/fig2_cell_vs_lung_cancer_%j.out
#SBATCH --time=06:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2 final comparison: cell vs. Xenium lung_cancer.
# Fully self-contained: generates the sim data (if not already present),
# QC-filters it, then runs all 4 count_distribution.py modes.
#
# Sim dataset: simulated single CELLS (no spatial aggregation), compared
# against Xenium lung cancer tissue (392-gene targeted panel).
# Current winning sim config (2026-08-25): log_mu=-2.5, theta=0.25,
# theta_jitter=0.15 -- see cell_vs_non_diseased_lung.sh for the full
# rationale (same sim config, different real reference).

# 2026-08-27: SIM_TAG now carries the per-modality batch_sigma finalized for
# the SHARED Figure 2 / Figure 3 dataset (cell 1.5, bin8 0.8, bin16 0.7,
# spot 0.3), tuned on the Figure 3 pre/post-Harmony demo then confirmed here
# to still match the real count distribution. count_distribution.py now
# defaults to --slice-id 5 and post-batch counts, so those flags are no
# longer passed per-call below. The matching clustering panel is
# code/clustering/<modality>_celltype_panel.sh (same SIM_TAG).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="log_mu_-2.5_theta_0.25_jitter0.15_bsigma15"
MODALITY="cell"
REAL_LABEL="lung_cancer"
SIM_RAW="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/figure_2/${MODALITY}_vs_${REAL_LABEL}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --base-gene-lognormal -2.5 0.7 \
        --theta 0.25 \
        --theta-jitter 0.15 \
        --batch-sigma 1.5 \
        --out-tag "${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
fi

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/full_panel"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/hvg_matched"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/qc_filtered"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

echo "[done] All 4 modes written under ${OUT_ROOT}"
