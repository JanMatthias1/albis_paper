#!/bin/bash
#SBATCH --job-name=fig2_larger_bin_vs_breast_cancer_visium_hd
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/larger_sphere/logs/fig2_larger_bin_vs_breast_cancer_visium_hd_%j.out
#SBATCH --time=24:00:00
#SBATCH --mem=400G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# larger_sphere counterpart of smaller_sphere/bin_vs_breast_cancer_visium_hd.sh
# (8um bins) -- grows the disc to cell's native sphere_r_um=6000 instead of
# the shared 2050um disc, holding 3D packing fraction fixed via
# n_cells = 600000 * (6000/2050)^3 ~ 15,043,310 (same target fraction as the
# smaller_sphere packing_pf0p04 config, ~4%).
#
# NEW territory as of 2026-09-17: unlike cell (already had a native-r6000
# baseline) and bin16um/spot (attempted the night before at 64G/8h and
# failed -- bin16um's matching n_cells task TIMED OUT, spot's OOM'd and
# TIMED OUT, see figure.md 2026-09-16), bin8um was never tried at r=6000 at
# all. Sized generously off that failure data: bin16um's 5.6M-cell task hit
# ~97GB RSS in 1h43m and its 15M-cell task didn't finish generating within
# 8h, so 400G/24h is a deliberately large first attempt, not a tuned
# estimate -- expect to revisit if it's still not enough.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/larger_sphere/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_log_mu_0.0_bsigma08_r6000"
MODALITY="bin"
REAL_LABEL="breast_cancer_visium_hd"
SIM_RAW="sim_paper/data/figure_2/larger_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/larger_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="sim_paper/data/figure_2/larger_sphere/plots/${MODALITY}_vs_${REAL_LABEL}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/larger_sphere/data/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 6000 \
            --n-cells 15043310 \
            --base-gene-lognormal 0.0 0.7 \
            --batch-sigma 0.8 \
            --sync-unaligned-seed \
            --out-tag "${SIM_TAG}"
    fi
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
            --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
    fi
    mkdir -p sim_paper/data/figure_2/larger_sphere/data
    mv "${NOISY_DIR}" "sim_paper/data/figure_2/larger_sphere/data/${SIM_TAG}"
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
