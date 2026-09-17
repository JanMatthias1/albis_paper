#!/bin/bash
#SBATCH --job-name=tune_banksy_cell
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/tune_banksy_cell_%A_%a.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#SBATCH --array=0-2
#
# One-off BANKSY parameter sweep for CELL's domain panel, mirroring the same
# sweep already run for bin (see figure.md, 2026-08-25). Default BANKSY
# settings (lambda=0.8, k_geom=15) collapsed cell's domain_true AND
# cell_type_true to exactly 0.0/0.0 on the strong-domain-mix dataset --
# testing whether that's a genuine limit or just under-tuned, same as the
# premature "BANKSY doesn't work for bin" conclusion that turned out to
# need more parameter search.
# 3 array tasks: 0=nbr_weight_decay=reciprocal, 1=max_m=0, 2=k_geom=30.
# Not a permanent script -- delete once cell_domain_panel.sh is updated
# with whichever config wins (or confirmed none do).

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs

# Two envs, matching 01_build_banksy_matrix.py's cross-env handoff design:
# the BANKSY build step needs the isolated sim-app-banksy stack; everything
# downstream (ari_vs_ground_truth.py) runs in the normal albis-tutorial
# env, same as misc/run_build_banksy_matrix.sh + the other panel scripts.
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"

cd /dcs04/hicks/data/Jan/sim_project

TAG="log_mu_-2.5_theta_0.25_strongdomainmix"
MODALITY="cell"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
CLUSTER_ROOT="sim_paper/data/clustering_${TAG}/${MODALITY}"

case "${SLURM_ARRAY_TASK_ID}" in
  0) VARIANT="reciprocal"; EXTRA_ARGS="--nbr-weight-decay reciprocal" ;;
  1) VARIANT="maxm0"; EXTRA_ARGS="--max-m 0" ;;
  2) VARIANT="kgeom30"; EXTRA_ARGS="--k-geom 30" ;;
esac

OUT_DIR="${CLUSTER_ROOT}/banksy_pca_harmony_qc_${VARIANT}"
echo "[banksy] variant=${VARIANT} args=${EXTRA_ARGS}"
"${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
    --modality "${MODALITY}" --packing-tag "${TAG}" --input "${SIM_QC}" \
    ${EXTRA_ARGS} \
    --output-dir "${OUT_DIR}"

echo "[ari] resolution-matched ARI recovery"
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/ari_vs_ground_truth.py \
    --modality "${MODALITY}" \
    --packing-tag "${TAG}" \
    --input "${OUT_DIR}/simulation_${MODALITY}_z_banksy_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/banksy_ari_recovery_${VARIANT}"

echo "[done] ${VARIANT} -> ${CLUSTER_ROOT}/banksy_ari_recovery_${VARIANT}/ari_summary_${MODALITY}.json"
