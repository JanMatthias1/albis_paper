#!/bin/bash
#SBATCH --job-name=fig2_bin8um_same_tissue_as_16um
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/smaller_sphere/logs/fig2_bin8um_same_tissue_as_16um_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# CANONICAL Figure 2 bin8um comparison as of 2026-09-23 (started as an
# exploratory check the same day; promoted by user decision because it matches
# real 8um Visium HD better on the panel-matched QC comparison: median
# counts/genes/zero-frac 47/21/0.953 vs breast 35/25/0.951, where the old
# log_mu=0.0 config gave 512/145/0.727). Replaces bin_vs_{breast_cancer,
# human_pancreas}_visium_hd.sh, archived with their data/plots under
# code/misc/count_distribution/misc/superseded_bin8um_logmu0_20260923/ and
# data/figure_2/misc/archive/bin8um_logmu0_superseded_20260923/.
#
# Question: could bin8um and bin16um share ONE underlying tissue? This takes
# bin16um's exact generate command (bin16um_vs_*_visium_hd.sh: log_mu=-2.5,
# theta 2.0, jitter 0.6, batch_sigma 0.7, same default seed) and changes only
# --bin-size-um 16 -> 8. The simulator bins the same full molecule stream
# (aggregate_molecules_to_grid_bins_2d_slices_window) and draws batch factors
# from seed-derived RNGs independent of bin size, so these 8um bins are an
# exact 2x2 subdivision of the canonical bin16um dataset's bins.
#
# The superseded 8um config used log_mu=0.0 / batch_sigma 0.8, i.e. a ~12x
# denser molecule tissue tuned separately from 16um.
# Compared against both real 8um Visium HD datasets; QC panels only (see
# realwindow note in the bin16um scripts -- pre-QC panels are ~77% empty).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/smaller_sphere/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bin8um_from16umcfg_log_mu_-2.5_jitter0.6_bsigma07"
MODALITY="bin"
SIM_RAW="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 2050 \
            --bin-size-um 8 \
            --base-gene-lognormal -2.5 0.7 \
            --theta 2.0 \
            --theta-jitter 0.6 \
            --batch-sigma 0.7 \
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

for REAL_LABEL in breast_cancer_visium_hd human_pancreas_visium_hd; do
    REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
    OUT_ROOT="sim_paper/data/figure_2/smaller_sphere/plots/bin8um_same_tissue_as_16um_vs_${REAL_LABEL}"

    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
        --output-dir "${OUT_ROOT}/qc_filtered"

    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
        --match-panel-size \
        --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

    echo "[done] ${OUT_ROOT}"
done
