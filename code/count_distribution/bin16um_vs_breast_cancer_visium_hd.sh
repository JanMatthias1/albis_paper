#!/bin/bash
#SBATCH --job-name=fig2_bin16um_vs_breast_cancer_visium_hd
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/fig2_bin16um_vs_breast_cancer_visium_hd_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2 final comparison: bin vs. Visium HD breast_cancer_visium_hd, at the
# 16um bin resolution (vs. the manuscript-baseline 8um used by
# bin_vs_breast_cancer_visium_hd.sh). Fully self-contained: generates the sim
# data (if not already present), QC-filters it, then runs all 4
# count_distribution.py modes.
#
# Sim dataset: simulated 16um BINS (real Visium HD ships 2um/8um/16um
# together, so this is a legitimate alternate real configuration, not just a
# tuning knob -- see generate_simulation_noisy.py --bin-size-um help and
# figure.md 2026-08-25 "bin's near-binary empty/non-empty bin sampling").
# Compared against real Visium HD breast cancer resampled to 16um
# (code/real_data_qc/misc/run_visium_hd_qc_16um.sh).
#
# Current winning sim config (2026-08-26, bin16um_logmu_sweep -- see
# figure.md / DATA_VERSIONS.md): packing_pf0p04 (sphere_r_um=2050, ~4% 3D
# packing, same as 8um -- packing is a 3D property independent of bin size),
# log_mu=-2.5 (deliberately more negative than 8um's winning log_mu=0.0,
# since a 16um bin covers ~4x the area of an 8um bin at fixed molecular
# density). theta/theta_jitter/noise_scale/marker_foldchange left at
# generate_simulation_noisy.py defaults (2.0/1.0/1.3/3.5), matching 8um --
# the 16um sweep only varied log_mu (and, in later rounds, theta/jitter,
# neither of which improved on the default).

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bin16um_log_mu_-2.5"
MODALITY="bin"
REAL_LABEL="breast_cancer_visium_hd_16um"
SIM_RAW="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/breast_cancer_visium_hd_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/figure_2/${MODALITY}_vs_${REAL_LABEL}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um 2050 \
        --bin-size-um 16 \
        --base-gene-lognormal -2.5 0.7 \
        --out-tag "${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
fi

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" --slice-id 5 \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/full_panel"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" --slice-id 5 \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/hvg_matched"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" --slice-id 5 \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/qc_filtered"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" --slice-id 5 \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

echo "[done] All 4 modes written under ${OUT_ROOT}"
