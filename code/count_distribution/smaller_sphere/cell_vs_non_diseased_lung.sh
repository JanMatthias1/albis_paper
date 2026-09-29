#!/bin/bash
#SBATCH --job-name=fig2_cell_vs_non_diseased_lung
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/smaller_sphere/logs/fig2_cell_vs_non_diseased_lung_%j.out
#SBATCH --time=06:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2, cell: simulated single cells vs Xenium non-diseased human lung (Xenium, 392-gene panel).
# The same simulation is compared with cell_vs_lung_cancer.sh (lung cancer).
# Settings: r = 2050 µm, 24,207 cells (the same physical tissue size as the
# bin/spot simulations), --base-gene-lognormal -2.5 0.7, theta 0.40,
# theta-jitter 0.15, batch_sigma 1.5 (set on the Figure 3 Harmony demo).
# log_mu and theta were chosen by a joint grid search against both Xenium
# references; theta 0.40 offsets the extra dispersion the batch effect adds.
# Figure 3's cell panels use the same data. Cells lie only inside the tissue,
# so all four modes are meaningful.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/smaller_sphere/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="log_mu_-2.5_theta_0.40_jitter0.15_bsigma15"
MODALITY="cell"
REAL_LABEL="non_diseased_lung"
SIM_RAW="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="sim_paper/data/figure_2/smaller_sphere/plots/${MODALITY}_vs_${REAL_LABEL}"

# generate_simulation_noisy.py only writes under data/noisy/<out-tag>/, so
# (re)generate + QC there and then move the directory into data/figure_2/.
if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 2050 \
            --n-cells 24207 \
            --base-gene-lognormal -2.5 0.7 \
            --theta 0.40 \
            --theta-jitter 0.15 \
            --batch-sigma 1.5 \
            --sync-unaligned-seed \
            --out-tag "${SIM_TAG}"
    fi
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/clustering/step00_qc_filter.py \
            --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
    fi
    mkdir -p sim_paper/data/figure_2/smaller_sphere/data
    mv "${NOISY_DIR}" "sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running step00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/step00_qc_filter.py \
        --modality "${MODALITY}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi

# full_panel: mean_variance/mean_dropout dispersion fit -- full gene panel,
# no HVG-matching (HVG selection would bias the dispersion estimate itself).
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/full_panel"

# hvg_matched: sparsity_summary/empty-row-rate -- panel-matched but NOT
# QC-filtered (QC-filtering sim here would hide the zero-inflation this
# stat exists to report).
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/hvg_matched"

# qc_filtered: total_counts (secondary use) with QC'd sim, full gene panel.
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/qc_filtered"

# qc_and_hvg_matched: total_counts/genes_per_cell -- both corrections
# together, the fair like-for-like for these two panels specifically.
"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

echo "[done] All 4 modes written under ${OUT_ROOT}"
