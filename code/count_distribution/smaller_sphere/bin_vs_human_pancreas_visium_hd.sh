#!/bin/bash
#SBATCH --job-name=fig2_bin_vs_human_pancreas_visium_hd
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/smaller_sphere/logs/fig2_bin_vs_human_pancreas_visium_hd_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2 final comparison: bin vs. Visium HD human_pancreas_visium_hd.
# Fully self-contained: generates the sim data (if not already present),
# QC-filters it, then runs all 4 count_distribution.py modes.
#
# Sim dataset: simulated 8um BINS (Visium-HD-like spatial aggregation),
# compared against real Visium HD human pancreas (whole-transcriptome panel).
# Current best sim config (2026-08-25, still under active tuning -- see
# bin_vs_breast_cancer_visium_hd.sh for the same config's rationale).

# 2026-08-27: SIM_TAG now carries the per-modality batch_sigma finalized for
# the SHARED Figure 2 / Figure 3 dataset (cell 1.5, bin8 0.8, bin16 0.7,
# spot 0.3), tuned on the Figure 3 pre/post-Harmony demo then confirmed here
# to still match the real count distribution. count_distribution.py now
# defaults to --slice-id 5 and post-batch counts, so those flags are no
# longer passed per-call below. The matching clustering panel is
# code/clustering/<modality>_celltype_panel.sh (same SIM_TAG).
#
# 2026-08-31: realwindow is now the default (the old sphere-scaled tight
# capture window is retired). generate_simulation_noisy.py no longer scales
# the window with --sphere-r-um; with --capture-window-um unset it uses the
# real instrument window (6.5 x 6.5 mm for Visium / Visium HD). The ~2050 um
# tissue disc sits inside that window with a wide empty border: off-tissue
# bins get domain_true / cell_type_true = "unassigned", obs["is_empty"] = True,
# and are dropped by 00_qc_filter.py (~77% of rows). Only the qc_filtered /
# qc_and_hvg_matched panels are meaningful for the figure; full_panel /
# hvg_matched are now empty-swamped diagnostics. Sim data is (re)generated
# under data/noisy/<tag>/ then moved into data/figure_2/<tag>/; pre-realwindow
# data is archived at data/figure_2_oldwindow_20260831/.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/smaller_sphere/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_log_mu_0.0_bsigma08"
MODALITY="bin"
REAL_LABEL="human_pancreas_visium_hd"
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
            --base-gene-lognormal 0.0 0.7 \
            --batch-sigma 0.8 \
            --sync-unaligned-seed \
            --out-tag "${SIM_TAG}"
    fi
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
            --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
    fi
    mkdir -p sim_paper/data/figure_2/smaller_sphere/data
    mv "${NOISY_DIR}" "sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
        --modality "${MODALITY}" --input "${SIM_RAW}" --output "${SIM_QC}"
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
