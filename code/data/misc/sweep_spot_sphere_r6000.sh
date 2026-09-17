#!/bin/bash
#SBATCH --job-name=sweep_spot_sphere_r6000
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_spot_sphere_r6000_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# spot counterpart of sweep_bin16um_sphere_r6000.sh -- see that script for
# the full rationale (growing the shared disc to cell's native sphere_r_um=
# 6000 instead of shrinking cell down to 2050, while holding 3D packing
# fraction fixed via n_cells ~ sphere_r_um^3 so bin/spot's already-solved
# empty-bin fix should carry over unchanged).
#
# n_cells grid, same 3 fractions as bin16um's sweep, at sphere_r_um=6000:
#   0.8%  -> ~3,008,662
#   1.5%  -> ~5,641,241
#   4.0%  -> ~15,043,310  (exact density match to current packing_pf0p04
#             spot config, n_cells=600000 @ sphere_r_um=2050)
#
# Dispersion knobs held at spot's current canonical config (log_mu=-2.0,
# theta=0.25, theta_jitter=0.10, batch_sigma=0.3 -- see
# spot_vs_{lymph_node,tonsil}_visium.sh). Compared against BOTH probe refs.
#
# Does NOT touch bin16um or cell -- see sweep_bin16um_sphere_r6000.sh /
# gen_cell_native_r6000.sh. If this lands well, promoting it means changing
# spot_vs_{lymph_node,tonsil}_visium.sh's --sphere-r-um/--n-cells (and
# spot_celltype_panel.sh / 3D_stair.py's spot entry, which share this tag)
# and regenerating -- not done automatically here.
#
# Tabulate with: sim_paper/code/data/misc/summary_spot_sphere_r6000.sh
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

N_CELLS_GRID=(3008662 5641241 15043310)
N_CELLS="${N_CELLS_GRID[${SLURM_ARRAY_TASK_ID:-0}]}"
LOG_MU=-2.0
THETA=0.25
JITTER=0.10
BATCH_SIGMA=0.3
SPHERE_R_UM=6000
MODALITY="spot"
TAG="spot_sphere_r6000_ncells${N_CELLS}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/spot_sphere_r6000"

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

for real in tonsil_visium lymph_node_visium; do
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
