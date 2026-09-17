#!/bin/bash
#SBATCH --job-name=gen_cell_native_r6000
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/gen_cell_native_r6000_%j.out
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Companion to sweep_bin16um_sphere_r6000.sh / sweep_spot_sphere_r6000.sh:
# regenerates cell at its ORIGINAL native geometry (sphere_r_um=6000,
# n_cells=600000 -- pre-2026-09-15 shrink) with cell's current dispersion
# knobs (log_mu=-2.3, theta=0.40, jitter=0.15, batch_sigma=1.5), so that if
# the bin16um/spot sphere-scale-up sweeps land well, there's a complete
# "everything at sphere_r_um=6000" picture to compare against tomorrow
# morning alongside sweep_cell_packing_r2050.sh's "everything at 2050"
# alternative.
#
# Not a sweep -- single point, no array. Full 4-mode count_distribution
# comparison against both Xenium refs (matching cell_vs_lung_cancer.sh /
# cell_vs_non_diseased_lung.sh's own structure) so the output is directly
# comparable to the current production cell_vs_{lung_cancer,
# non_diseased_lung} panels.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

LOG_MU=-2.3
THETA=0.40
JITTER=0.15
BATCH_SIGMA=1.5
SPHERE_R_UM=6000
N_CELLS=600000
MODALITY="cell"
TAG="cell_native_r6000"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/cell_native_r6000"

echo "[config] n_cells=${N_CELLS} sphere_r_um=${SPHERE_R_UM} log_mu=${LOG_MU} theta=${THETA} jitter=${JITTER} batch_sigma=${BATCH_SIGMA} tag=${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um "${SPHERE_R_UM}" \
        --n-cells "${N_CELLS}" \
        --base-gene-lognormal "${LOG_MU}" 0.7 \
        --theta "${THETA}" \
        --theta-jitter "${JITTER}" \
        --batch-sigma "${BATCH_SIGMA}" \
        --sync-unaligned-seed \
        --out-tag "${TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
fi

for real in lung_cancer non_diseased_lung; do
    REAL_INPUT="sim_paper/data/real_data_qc/${real}/${real}_qc.h5ad"

    echo "[compare] ${TAG} vs ${real} (qc_filtered)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "${REAL_INPUT}" --compare-label "${real}" \
        --output-dir "${OUT_ROOT}/${TAG}_qc_filtered_vs_${real}"

    echo "[compare] ${TAG} vs ${real} (hvg_matched)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "${REAL_INPUT}" --compare-label "${real}" \
        --match-panel-size \
        --output-dir "${OUT_ROOT}/${TAG}_hvg_matched_vs_${real}"
done

echo "[done] ${TAG}"
