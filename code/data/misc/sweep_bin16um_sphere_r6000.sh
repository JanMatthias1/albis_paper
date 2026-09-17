#!/bin/bash
#SBATCH --job-name=sweep_bin16um_sphere_r6000
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_bin16um_sphere_r6000_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Alternative fix for the cell-vs-xenium regression (see
# sweep_cell_packing_r2050.sh for the other approach): instead of shrinking
# cell down to bin16um/spot's small disc, grow bin16um's disc UP to cell's
# ORIGINAL native scale (sphere_r_um=6000) so cell can just revert to its
# proven-good config untouched (gen_cell_native_r6000.sh) while bin16um/spot
# meet it there instead.
#
# WHY THIS SHOULD BE SAFE: bin16um's --n-cells/--sphere-r-um were tuned to
# packing_pf0p04 (n_cells=600000 @ sphere_r_um=2050) specifically to fix a
# 66.3%-empty-bins problem at the old 6000um/600k geometry (see
# generate_simulation_noisy.py header comment, ~line 118). That fix is about
# 3D PACKING FRACTION (n_cells / sphere volume), not absolute sphere size --
# cell_r_mean (molecule-spillover radius) is a fixed absolute default that
# does NOT scale with --sphere-r-um, and n_cells scaling as sphere_r_um^3
# holds the mean inter-cell spacing constant regardless of R. So the same
# fraction realized at R=6000 with proportionally more cells should reproduce
# the same empty-bin/overshoot behaviour as at R=2050 -- this sweep is the
# check on that claim, not a re-tune from scratch.
#
# n_cells grid targets 3 fractions spanning the previously-documented
# "best compromise" range (0.8-1.5% empty-bins-vs-overshoot, generate_
# simulation_noisy.py header) up to the exact match of the CURRENT
# packing_pf0p04 config (~4%), all at sphere_r_um=6000:
#   0.8%  -> n_cells = 600000 * (6000/2050)^3 * (0.008/0.04) ~ 3,008,662
#   1.5%  -> ~5,641,241
#   4.0%  -> ~15,043,310  (exact density match to current bin16um config)
# Cost check: an existing 600k-cell run at this same code path took ~5 min /
# 2.3GB RAM (job 35718698) -- even a generous non-linear safety margin should
# comfortably fit the 64G/8h requested here for up to ~15M cells.
#
# Dispersion knobs held at bin16um's current canonical config (log_mu=-2.5,
# batch_sigma=0.7, bin_size_um=16 -- see bin16um_vs_breast_cancer_visium_hd.sh
# / bin16um_vs_human_pancreas_visium_hd.sh). Compared against BOTH real refs.
#
# Does NOT touch spot or cell -- see sweep_spot_sphere_r6000.sh /
# gen_cell_native_r6000.sh for their counterparts. If this lands well,
# promoting it means changing bin16um_vs_*.sh's --sphere-r-um/--n-cells and
# regenerating the shared Fig2/Fig3 bin16um dataset -- not done automatically
# here.
#
# Tabulate with: sim_paper/code/data/misc/summary_bin16um_sphere_r6000.sh
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

N_CELLS_GRID=(3008662 5641241 15043310)
N_CELLS="${N_CELLS_GRID[${SLURM_ARRAY_TASK_ID:-0}]}"
LOG_MU=-2.5
BATCH_SIGMA=0.7
SPHERE_R_UM=6000
BIN_SIZE_UM=16
MODALITY="bin"
TAG="bin16um_sphere_r6000_ncells${N_CELLS}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/bin16um_sphere_r6000"

echo "[config] n_cells=${N_CELLS} sphere_r_um=${SPHERE_R_UM} bin_size_um=${BIN_SIZE_UM} log_mu=${LOG_MU} batch_sigma=${BATCH_SIGMA} tag=${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um "${SPHERE_R_UM}" \
        --n-cells "${N_CELLS}" \
        --bin-size-um "${BIN_SIZE_UM}" \
        --base-gene-lognormal "${LOG_MU}" 0.7 \
        --batch-sigma "${BATCH_SIGMA}" \
        --sync-unaligned-seed \
        --out-tag "${TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
fi

for real in breast_cancer_visium_hd human_pancreas_visium_hd; do
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
