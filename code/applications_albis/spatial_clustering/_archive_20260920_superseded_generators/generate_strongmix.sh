#!/bin/bash
# ARCHIVED 2026-09-20: superseded by code/clustering/strong_mix/generate_strong_mix_{cell,bin16um,spot}.sh
# + stagate_inputs.py's DATASETS dict (3D_stagate.py's --dataset choices no longer
# include this script's tags, so it's unreachable). Kept for history only --
# its cell config (log_mu=-2.3, default sphere_r_um=6000/n_cells=600000) is
# STALE, predating the 2026-09-17 dispersion retune and diameter match to
# Figure 2 (log_mu=-2.5, sphere_r_um=2050, n_cells=24207). Do not rerun.
#SBATCH --job-name=stagate_strongmix
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering/logs_stagate_3D/strongmix_%A_%a.out
#SBATCH --array=0-5
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 4B (STAGATE) data: the strong-domain-mix half of the 2x2
# {weak|strong domain mix} x {tuned|very-low batch} grid. The weak half is
# the Figure 2 tags (family 1, already on disk) + generate_figure2_lowbatch.sh
# (family 3). This script builds:
#   family 2  (strong mix, TUNED batch):   tasks 0-2   bsigma 0.7 / 0.3 / 1.5
#   family 4  (strong mix, VERY LOW batch): tasks 3-5   bsigma 0.05
#
# Every config is IDENTICAL to the matching Figure 2 tag except (a) the
# --strong-domain-mix flag and (b) --batch-sigma, so the only things varying
# across the whole 2x2 are the domain-mix flag and the batch strength.
# Current realwindow default (no window scaling) + current post-batch QC
# (00_qc_filter.py), so it is consistent with families 1 and 3.
#
#   bin -> bin16um: log_mu=-2.5, theta=2.0, jitter=1.0, sphere_r_um=2050, bin_size_um=16
#   spot:           log_mu=-2.5, theta=2.0, jitter=1.0, sphere_r_um=2050
#   cell:           log_mu=-2.3, theta=0.40, jitter=0.15, sphere_r_um=default 6000
#
# generate_simulation_noisy.py only writes to data/noisy/<out-tag>/, so this
# generates+QCs there and then moves the dir into
# data/figure_4/spatial_clustering/sim_data/<tag>/ (all Figure-4B-only data
# lives there). Safe to rerun: skips any tag already present at the target.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering/logs_stagate_3D
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

#            fam2 (tuned batch)                fam4 (very low batch)
MODALITIES=(bin        spot       cell         bin        spot       cell)
LOG_MUS=(   -2.5       -2.5       -2.3         -2.5       -2.5       -2.3)
BIN_SIZES=( 16         -          -            16         -          -)
THETAS=(    2.0        2.0        0.40         2.0        2.0        0.40)
JITTERS=(   1.0        1.0        0.15         1.0        1.0        0.15)
BATCH_SIGMAS=(0.7      0.3        1.5          0.05       0.05       0.05)
TAGS=(
  packing_pf0p04_bin16um_log_mu_-2.5_bsigma07_strongmix
  packing_pf0p04_log_mu_-2.5_bsigma03_strongmix
  log_mu_-2.3_theta_0.40_jitter0.15_bsigma15_strongmix
  packing_pf0p04_bin16um_log_mu_-2.5_bsigma005_strongmix
  packing_pf0p04_log_mu_-2.5_bsigma005_strongmix
  log_mu_-2.3_theta_0.40_jitter0.15_bsigma005_strongmix
)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
MODALITY="${MODALITIES[$IDX]}"
LOG_MU="${LOG_MUS[$IDX]}"
BIN_SIZE="${BIN_SIZES[$IDX]}"
THETA="${THETAS[$IDX]}"
JITTER="${JITTERS[$IDX]}"
BATCH_SIGMA="${BATCH_SIGMAS[$IDX]}"
TAG="${TAGS[$IDX]}"

NOISY_DIR="sim_paper/data/noisy/${TAG}"
DEST_DIR="sim_paper/data/figure_4/spatial_clustering/sim_data/${TAG}"

echo "[config] idx=${IDX} modality=${MODALITY} log_mu=${LOG_MU} bin_size=${BIN_SIZE} theta=${THETA} jitter=${JITTER} batch_sigma=${BATCH_SIGMA} strong_domain_mix=True tag=${TAG}"

if [[ -f "${DEST_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
    echo "[skip] ${DEST_DIR} already complete"
    exit 0
fi

GEN_ARGS=(--modality "${MODALITY}" --base-gene-lognormal "${LOG_MU}" 0.7 --batch-sigma "${BATCH_SIGMA}" --strong-domain-mix --out-tag "${TAG}")
if [[ "${MODALITY}" == "bin" || "${MODALITY}" == "spot" ]]; then
    GEN_ARGS+=(--sphere-r-um 2050)
fi
if [[ "${MODALITY}" == "cell" ]]; then
    GEN_ARGS+=(--theta "${THETA}" --theta-jitter "${JITTER}")
fi
if [[ "${BIN_SIZE}" != "-" ]]; then
    GEN_ARGS+=(--bin-size-um "${BIN_SIZE}")
fi

if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
    echo "[generate] ${TAG}"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py "${GEN_ARGS[@]}"
fi

if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
    echo "[qc] ${TAG}"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
fi

echo "[move] ${NOISY_DIR} -> ${DEST_DIR}"
mkdir -p sim_paper/data/figure_4/spatial_clustering/sim_data
mv "${NOISY_DIR}" "${DEST_DIR}"

echo "[done] ${TAG}"
