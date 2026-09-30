#!/bin/bash
#SBATCH --job-name=fig2_lowbatch
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering/logs_stagate_3D/fig2_lowbatch_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 4B (STAGATE), combination #3a: the 3 canonical Figure 2 configs
# regenerated with a VERY LOW batch effect (batch_sigma = 0.05) and nothing
# else changed. Pairs against the tuned-batch Figure 2 tags (cell 1.5 /
# bin16 0.7 / spot 0.3) to isolate what the batch effect does to STAGATE's
# 3D-vs-2D domain recovery on weak-domain-mix data.
#
# Params copied verbatim from code/misc/data/misc/generate_figure2_batch_tuned.sh
# (which read them off each existing figure_2/<tag>/simulation_*.h5ad
# uns['sim_params']):
#   bin 16um: log_mu=-2.5, theta=2.0, jitter=1.0, sphere_r_um=2050, bin_size_um=16
#   spot:     log_mu=-2.5, theta=2.0, jitter=1.0, sphere_r_um=2050
#   cell:     log_mu=-2.3, theta=0.40, jitter=0.15, sphere_r_um=default 6000
# Only --batch-sigma differs (-> 0.05). Uses the current (post-batch) QC in
# step00_qc_filter.py and the current realwindow default in
# generate_simulation_noisy.py, so it is consistent with the family-1 tags.
#
# generate_simulation_noisy.py writes to data/noisy/<out-tag>/ only; this
# script generates+QCs there then moves the dir into
# data/figure_4/spatial_clustering/sim_data/<out-tag>/ (all Figure-4B-only data
# lives there). Safe to rerun: skips any tag already present at the target.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering/logs_stagate_3D
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

MODALITIES=(bin                                        spot                              cell)
LOG_MUS=(-2.5                                          -2.5                              -2.3)
BIN_SIZES=(16                                          -                                 -)
THETAS=(2.0                                            2.0                               0.40)
JITTERS=(1.0                                           1.0                               0.15)
BATCH_SIGMA=0.05
TAGS=(
  packing_pf0p04_bin16um_log_mu_-2.5_bsigma005
  packing_pf0p04_log_mu_-2.5_bsigma005
  log_mu_-2.3_theta_0.40_jitter0.15_bsigma005
)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
MODALITY="${MODALITIES[$IDX]}"
LOG_MU="${LOG_MUS[$IDX]}"
BIN_SIZE="${BIN_SIZES[$IDX]}"
THETA="${THETAS[$IDX]}"
JITTER="${JITTERS[$IDX]}"
TAG="${TAGS[$IDX]}"

NOISY_DIR="sim_paper/data/noisy/${TAG}"
DEST_DIR="sim_paper/data/figure_4/spatial_clustering/sim_data/${TAG}"

echo "[config] idx=${IDX} modality=${MODALITY} log_mu=${LOG_MU} bin_size=${BIN_SIZE} theta=${THETA} jitter=${JITTER} batch_sigma=${BATCH_SIGMA} tag=${TAG}"

if [[ -f "${DEST_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
    echo "[skip] ${DEST_DIR} already complete"
    exit 0
fi

GEN_ARGS=(--modality "${MODALITY}" --base-gene-lognormal "${LOG_MU}" 0.7 --batch-sigma "${BATCH_SIGMA}" --out-tag "${TAG}")
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
    "${PYTHON_BIN}" sim_paper/code/clustering/step00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
fi

echo "[move] ${NOISY_DIR} -> ${DEST_DIR}"
mkdir -p sim_paper/data/figure_4/spatial_clustering/sim_data
mv "${NOISY_DIR}" "${DEST_DIR}"

echo "[done] ${TAG}"
