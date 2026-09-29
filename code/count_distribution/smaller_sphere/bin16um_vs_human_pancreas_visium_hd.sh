#!/bin/bash
#SBATCH --job-name=fig2_bin16um_vs_human_pancreas_visium_hd
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/smaller_sphere/logs/fig2_bin16um_vs_human_pancreas_visium_hd_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2, bin 16 µm: simulated 16 µm bins vs Visium HD human pancreas, resampled to 16 µm bins
# (Visium HD ships 2/8/16 µm; QC by real_data_qc/). The same simulation is
# compared with bin16um_vs_breast_cancer_visium_hd.sh (breast cancer).
# Settings: r = 2050 µm, --base-gene-lognormal -2.5 0.7, theta 2.0,
# theta-jitter 0.6 (larger jitter floors some genes to extreme dispersion),
# batch_sigma 0.7 (set on the Figure 3 Harmony demo). Figure 3's bin16um panels
# use the same data.
# The simulation covers the full 6.5 × 6.5 mm window, so most bins are
# off-tissue before QC: use the qc_filtered / qc_and_hvg_matched outputs.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/smaller_sphere/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma07"
MODALITY="bin"
REAL_LABEL="human_pancreas_visium_hd_16um"
SIM_RAW="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/human_pancreas_visium_hd_qc.h5ad"
OUT_ROOT="sim_paper/data/figure_2/smaller_sphere/plots/${MODALITY}_vs_${REAL_LABEL}"

# generate_simulation_noisy.py only writes under data/noisy/<out-tag>/, so
# (re)generate + QC there and then move the directory into data/figure_2/.
if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 2050 \
            --bin-size-um 16 \
            --base-gene-lognormal -2.5 0.7 \
            --theta 2.0 \
            --theta-jitter 0.6 \
            --batch-sigma 0.7 \
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
