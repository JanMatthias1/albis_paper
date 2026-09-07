#!/bin/bash
#SBATCH --job-name=tune_banksy_bin_kgeom
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/tune_banksy_bin_kgeom_%A_%a.out
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#SBATCH --array=0-2
#
# Follow-up to the 2026-08-25 bin BANKSY sweep (default/reciprocal/max_m=0/
# k_geom=30, see figure.md), which found k_geom=30 as the best domain_true
# ARI so far (0.093 vs. plain's 0.048, both against
# packing_pf0p04_bsigma05_strongdomainmix) -- but the resolution binary
# search for k_geom=30 plateaued at exactly 10 clusters all the way down to
# the res_lo=0.01 floor (flat across the last 9/15 trials), never reaching
# the target k=6. That's a graph-fragmentation floor, not a resolution-
# range problem (widening res_lo below 0.01 won't unstick a value that's
# already flat), so the fix tested here is more spatial neighbors (higher
# k_geom) to bridge the residual small clusters into the main structure --
# same mechanism as the earlier QC-filter fix for bin's original 189-cluster
# fragmentation.
#
# 3 array tasks: 0=k_geom=50, 1=k_geom=75, 2=k_geom=100.
# Not a permanent script -- delete once bin_domain_panel.sh is updated with
# whichever config wins (or confirmed none do / achieved_k still misses target).

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs

# Two envs, matching 01_build_banksy_matrix.py's cross-env handoff design:
# the BANKSY build step needs the isolated sim-app-banksy stack; everything
# downstream (ari_vs_ground_truth.py) runs in the normal albis-tutorial
# env, same as the other panel/tuning scripts.
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"

cd /dcs04/hicks/data/Jan/sim_project

TAG="packing_pf0p04_bsigma05_strongdomainmix"
MODALITY="bin"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
CLUSTER_ROOT="sim_paper/data/clustering_${TAG}/${MODALITY}"

case "${SLURM_ARRAY_TASK_ID}" in
  0) VARIANT="kgeom50";  EXTRA_ARGS="--k-geom 50" ;;
  1) VARIANT="kgeom75";  EXTRA_ARGS="--k-geom 75" ;;
  2) VARIANT="kgeom100"; EXTRA_ARGS="--k-geom 100" ;;
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
