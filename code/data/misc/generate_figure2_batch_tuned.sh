#!/bin/bash
#SBATCH --job-name=fig2_batch_tuned
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/fig2_batch_tuned_%A_%a.out
#SBATCH --array=0-3
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Regenerates each modality's already-decided Figure 2 winning config,
# swapping in the per-modality batch_sigma finalized 2026-08-27 for the
# SHARED Figure 2 / Figure 3 dataset: cell 1.5, bin 8um 0.8, bin 16um 0.7,
# spot 0.3. These were picked on Figure 3's pre/post-Harmony slice-separation
# demo (batch_sigma_sweep.sh + batch_sigma_opt_round2.sh) and then confirmed
# against the real count distributions here -- superseding the earlier
# flat-0.22 and the 0.22/0.5/0.5/0.15 "batch_tuned" values.
# Everything else (theta/theta_jitter/base_gene_lognormal/sphere_r_um/
# bin_size_um) is held identical to the existing winning tag, confirmed by
# reading uns['sim_params'] off each existing figure_2/<tag>/simulation_*.h5ad:
#   bin 8um:  packing_pf0p04_log_mu_0.0                (theta=2.0, jitter=1.0, log_mu=0.0,  sphere_r_um=2050)
#   bin 16um: packing_pf0p04_bin16um_log_mu_-2.5        (theta=2.0, jitter=1.0, log_mu=-2.5, sphere_r_um=2050, bin_size_um=16)
#   spot:     packing_pf0p04_log_mu_-2.5                (theta=2.0, jitter=1.0, log_mu=-2.5, sphere_r_um=2050)
#   cell:     log_mu_-2.3_theta_0.40_jitter0.15         (theta=0.40, jitter=0.15, log_mu=-2.3, sphere_r_um=default 6000 -- not packing-affected; theta 0.40 pre-compensates the batch-effect dispersion drop, see figure.md 2026-08-27)
#
# generate_simulation_noisy.py writes to data/noisy/<out-tag>/ only (its
# BASE_DATA_DIR is hardcoded) -- this script generates+QCs there, then moves
# the whole directory into data/figure_2/<out-tag>/ itself, so this never
# hits the figure_2-path generate/QC mismatch documented in figure.md
# 2026-08-26 (the REAL_INPUT filename bug was separate, but the same root
# cause: the vs_*.sh compare scripts' own generate/QC fallback branches
# write to data/noisy/ while checking data/figure_2/, so they'd silently
# fail on a genuinely-missing tag). Safe to rerun: skips generate/QC/move
# for any tag already present under data/figure_2/.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

MODALITIES=(bin        bin                    spot       cell)
LOG_MUS=(0.0            -2.5                   -2.5       -2.3)
BIN_SIZES=(8            16                     -          -)
THETAS=(2.0             2.0                    2.0        0.40)
JITTERS=(1.0            1.0                    1.0        0.15)
BATCH_SIGMAS=(0.8       0.7                    0.3        1.5)
TAGS=(
  packing_pf0p04_log_mu_0.0_bsigma08
  packing_pf0p04_bin16um_log_mu_-2.5_bsigma07
  packing_pf0p04_log_mu_-2.5_bsigma03
  log_mu_-2.3_theta_0.40_jitter0.15_bsigma15
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
FIG2_DIR="sim_paper/data/figure_2/${TAG}"

echo "[config] idx=${IDX} modality=${MODALITY} log_mu=${LOG_MU} bin_size=${BIN_SIZE} theta=${THETA} jitter=${JITTER} batch_sigma=${BATCH_SIGMA} tag=${TAG}"

if [[ -f "${FIG2_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
    echo "[skip] ${FIG2_DIR} already complete"
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
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
fi

echo "[move] ${NOISY_DIR} -> ${FIG2_DIR}"
mkdir -p sim_paper/data/figure_2
mv "${NOISY_DIR}" "${FIG2_DIR}"

echo "[done] ${TAG}"
