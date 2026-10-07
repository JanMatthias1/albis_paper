#!/bin/bash
#SBATCH --job-name=fig2_spot_vs_lymph_node_visium
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/smaller_sphere/logs/fig2_spot_vs_lymph_node_visium_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2, spot: simulated Visium-like spots (100 µm spacing, 27.5 µm capture
# radius) vs Visium human lymph node (CytAssist FFPE probe panel, 18k genes;
# QC by real_data_qc/). The same simulation is compared with the second probe
# reference in spot_vs_tonsil_visium.sh (tonsil); the spot settings were tuned against both
# (composite of |ln(sim/real)| over θ̂, median total counts and zero fraction).
# Settings: r = 2050 µm, --base-gene-lognormal -2.25 1.0, theta 0.25,
# theta-jitter 0.10 (larger jitter floors many genes to extreme dispersion),
# batch_sigma 0.3 (set on the Figure 3 Harmony demo). Figure 3's spot panels
# use the same data.
# The simulation covers the full 6.5 × 6.5 mm Visium window, so most spots are
# off-tissue before QC: use the qc_filtered / qc_and_hvg_matched outputs.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/smaller_sphere/logs
source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_bsigma03"
MODALITY="spot"
REAL_LABEL="lymph_node_visium"
SIM_RAW="albis_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="albis_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="albis_paper/data/noisy/${SIM_TAG}"
REAL_INPUT="albis_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="albis_paper/data/figure_2/smaller_sphere/plots/${MODALITY}_vs_${REAL_LABEL}"

# generate_simulation_noisy.py only writes under data/noisy/<out-tag>/, so
# (re)generate + QC there and then move the directory into data/figure_2/.
if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" albis_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 2050 \
            --base-gene-lognormal -2.25 1.0 \
            --theta 0.25 \
            --theta-jitter 0.10 \
            --batch-sigma 0.3 \
            --sync-unaligned-seed \
            --out-tag "${SIM_TAG}"
    fi
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
        "${PYTHON_BIN}" albis_paper/code/clustering/step00_qc_filter.py \
            --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
    fi
    mkdir -p albis_paper/data/figure_2/smaller_sphere/data
    mv "${NOISY_DIR}" "albis_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running step00_qc_filter.py"
    "${PYTHON_BIN}" albis_paper/code/clustering/step00_qc_filter.py \
        --modality "${MODALITY}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi

if [[ ! -f "${REAL_INPUT}" ]]; then
    echo "ERROR: ${REAL_INPUT} not found -- run real_data_qc/run_visium_lymph_qc.sh first" >&2
    exit 1
fi

"${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/full_panel"

"${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/hvg_matched"

"${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/qc_filtered"

"${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

echo "[done] All 4 modes written under ${OUT_ROOT}"
