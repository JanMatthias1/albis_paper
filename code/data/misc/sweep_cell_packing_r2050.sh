#!/bin/bash
#SBATCH --job-name=sweep_cell_packing_r2050
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_cell_packing_r2050_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# cell --n-cells sweep at the fixed sphere_r_um=2050 disc (cross_tech_stair's
# shared cell/bin16um/spot physical scale). Dispersion knobs held at cell's
# current canonical config (log_mu=-2.3, theta=0.40, jitter=0.15,
# batch_sigma=1.5).
#
# WHY: on 2026-09-15 cell's sphere was shrunk from 6000um/600000 cells to
# 2050um/24207 cells to match bin16um/spot's disc, scaling n_cells by
# (2050/6000)^3 ~ 0.04 to hold cell's ORIGINAL ~0.16% packing fraction fixed
# (non-linear -- see figure.md 2026-09-16). theta_hat/sparsity are basically
# unchanged (0.08/92.5% now vs 0.08/92.15% before), but the per-slice sample
# feeding count_distribution.py's comparison collapsed from ~60,000 cells to
# 3,544 (slice_id=5) -- a ~17x smaller sample for every Figure 2 QC plot,
# which is the likely source of the "cell vs xenium isn't as good as before"
# regression reported by the user.
#
# This is one of two parallel approaches being tried tonight (see
# sweep_bin16um_sphere_r6000.sh / sweep_spot_sphere_r6000.sh /
# gen_cell_native_r6000.sh for the other: reverting cell to its native
# 6000um/600k config and scaling bin/spot UP to match instead). This script
# keeps cell's disc small and just adds more cells; that one keeps
# everything's density fixed and grows the shared disc instead.
#
# Sweeps --n-cells UP at the fixed 2050um sphere (raising packing fraction
# above cell's original 0.16%, toward bin/spot's ~4% packing_pf0p04
# convention). Does NOT touch bin/spot.
#
# Caveat: raising packing fraction can cause neighbouring cells' molecule
# clouds to overlap (this is what caused bin's 5-9x count/bin overshoot at
# high packing, generate_simulation_noisy.py header comment) -- possible
# this shows up here too as an upper bound on how far --n-cells can go
# before cell's per-cell count distribution overshoots real Xenium. That's
# exactly what this sweep is meant to reveal.
#
# Two comparison modes per (n_cells, ref): qc_filtered (full 556-gene panel,
# QC'd -- what the figure reads cell's mean_variance/mean_dropout from) and
# hvg_matched (apples-to-apples theta_hat/zero_frac vs the real ref's panel
# size).
#
# Baseline (n_cells=24207, already generated under data/figure_2/) is NOT
# regenerated here -- summary_cell_packing_r2050.sh pulls it straight from
# data/count_distribution/figure_2/cell_vs_{lung_cancer,non_diseased_lung}/.
#
# Tabulate with: sim_paper/code/data/misc/summary_cell_packing_r2050.sh
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

N_CELLS_GRID=(75000 200000 600000)
N_CELLS="${N_CELLS_GRID[${SLURM_ARRAY_TASK_ID:-0}]}"
LOG_MU=-2.3
THETA=0.40
JITTER=0.15
BATCH_SIGMA=1.5
SPHERE_R_UM=2050
MODALITY="cell"
TAG="cell_packing_sweep_ncells${N_CELLS}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/cell_packing_sweep_r2050"

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
