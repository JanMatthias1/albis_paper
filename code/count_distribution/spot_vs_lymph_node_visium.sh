#!/bin/bash
#SBATCH --job-name=fig2_spot_vs_lymph_node_visium
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/fig2_spot_vs_lymph_node_visium_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2 second real-Visium comparison: spot vs. Visium V2 Human Lymph
# Node (probe-based CytAssist FFPE -- see real_data_qc/run_visium_lymph_qc.sh
# for provenance), run alongside spot_vs_breast_cancer_visium.sh
# (whole-transcriptome, fresh-frozen) so the "spot" count-distribution match
# is checked against two different capture chemistries, not just one tissue.
# Same SIM_TAG / sim data as spot_vs_breast_cancer_visium.sh -- only the real
# reference changes -- so this is fully self-contained (generates the sim
# data if not already present) but in practice reuses what that script
# already built under data/figure_2/.
#
# Prereq: real_data_qc/run_visium_lymph_qc.sh (writes
# data/real_data_qc/lymph_node_visium/lymph_node_visium_qc.h5ad).
#
# 2026-09-04: initial run, ad hoc (not yet via sbatch) -- see
# data/count_distribution/figure_2/spot_vs_lymph_node_visium/*/comparison_summary.json.
# qc_filtered (no panel match): total_counts median ratio sim/real = 0.05
# (sim ~20x lower than breast_cancer_visium's 6.7x gap -- lymph node's probe
# panel is far deeper per spot than WTA); theta_hat ratio 0.59 (real less
# dispersed than sim here, opposite of breast_cancer's near-1:1 match).
# qc_and_hvg_matched (top-556 HVGs): total_counts ratio improves to 0.59,
# genes_per_cell median sim 358 vs real 394 (close), but theta_hat ratio
# flips to 1.87 (sim now more dispersed than real). Open question for
# figure.md: whether tuned batch_sigma=0.3 (fit against breast_cancer_visium)
# still holds against a probe-based reference, or needs a separate check.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_log_mu_-2.5_bsigma03"
MODALITY="spot"
REAL_LABEL="lymph_node_visium"
SIM_RAW="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/figure_2/${MODALITY}_vs_${REAL_LABEL}"

# generate_simulation_noisy.py only writes under data/noisy/<out-tag>/, so
# (re)generate + QC there and then move the directory into data/figure_2/.
if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 2050 \
            --base-gene-lognormal -2.5 0.7 \
            --batch-sigma 0.3 \
            --out-tag "${SIM_TAG}"
    fi
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
            --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
    fi
    mkdir -p sim_paper/data/figure_2
    mv "${NOISY_DIR}" "sim_paper/data/figure_2/${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
        --modality "${MODALITY}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi

if [[ ! -f "${REAL_INPUT}" ]]; then
    echo "ERROR: ${REAL_INPUT} not found -- run real_data_qc/run_visium_lymph_qc.sh first" >&2
    exit 1
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
