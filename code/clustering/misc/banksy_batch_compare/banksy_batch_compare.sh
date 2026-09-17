#!/bin/bash
#SBATCH --job-name=banksy_batch_cmp
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/banksy_batch_cmp_%A_%a.out
#SBATCH --array=0-5%6
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# BANKSY domain recovery vs per-slice batch effect -- the BANKSY counterpart of
# the STAGATE Figure 4B panel (plot_stagate_fig4b.py). At a FIXED (lambda=0.5,
# k_geom=200 -- the winner of banksy_lambda_kgeom_sweep.sh on near-batch-free
# data), for each resolution, run BANKSY on:
#
#   no batch effect   --use-pre-batch  (X <- layers['counts_pre_batch'])
#   tuned batch        the strong-mix h5ad as-is, per-resolution tuned sigma
#                      (cell 1.5 / bin16um 0.7 / spot 0.3)
#
# WHY: the lambda/k_geom sweep's clean cell result (domain ARI 0.76) turned out
# to be on a stale config with batch_sigma 0.22; re-running on the current
# tuned-batch strong-mix tags gave domain ARI ~0.02 with slice-leakage ARI ~1.0
# (clusters == slice_id). This isolates the batch effect as the cause and lets
# us state, side by side with STAGATE, whether BANKSY's Harmony step buys any
# batch tolerance (it appears not to at the tuned sigma).
#
# Same strong-mix datasets Figure 4B STAGATE uses, so the two panels are on
# byte-identical data. spot's strong-mix h5ad is on the older log_mu_-2.5 count
# config -- acceptable here since this measures batch structure, not counts
# (same rationale as 3D_stagate.py's spot note).
#
# Per task: (1) 01_build_banksy_matrix.py [sim-app-banksy env], (2)
# 02_leiden_resolution_sweep.py + (3) composition_recovery.py --h5ad [albis-tutorial],
# --stagger-scale 5, --max-m 1, --skip-umap. composition_recovery.py reports
# ARI(clusters, slice_id) as the leakage flag.
#
# Outputs under data/figure_3/banksy_batch_compare/:
#   <mod>/<batch>/banksy_matrix/..._banksy_pca_harmony_qc.h5ad
#   <mod>/<batch>/ari/ari_summary_<mod>.json
#   scores/composition_recovery_<mod>_<batch>.json
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

SIM_DATA_DIR="sim_paper/data/figure_4/spatial_clustering/sim_data"
LAM=0.5
KG=200

# (modality, strong-mix tag)
MODS=(cell bin spot)
TAGS=(
  "log_mu_-2.3_theta_0.40_jitter0.15_bsigma15_strongmix"
  "packing_pf0p04_bin16um_log_mu_-2.5_bsigma07_strongmix"
  "packing_pf0p04_log_mu_-2.5_bsigma03_strongmix"
)
BATCHES=(prebatch tuned)

i="${SLURM_ARRAY_TASK_ID:-0}"
MOD="${MODS[$(( i / 2 ))]}"
TAG="${TAGS[$(( i / 2 ))]}"
BATCH="${BATCHES[$(( i % 2 ))]}"

SIM_QC="${SIM_DATA_DIR}/${TAG}/simulation_${MOD}_z_qc.h5ad"
SWEEP="sim_paper/data/figure_3/banksy_batch_compare"
RUN="${SWEEP}/${MOD}/${BATCH}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_${MOD}_z_banksy_pca_harmony_qc.h5ad"
ARI_H5AD="${RUN}/ari/simulation_${MOD}_z_ari_recovery.h5ad"

PREBATCH_FLAG=()
[[ "${BATCH}" == "prebatch" ]] && PREBATCH_FLAG=(--use-pre-batch)

echo "[task ${i}] MOD=${MOD} batch=${BATCH} lambda=${LAM} k_geom=${KG}"
echo "           input=${SIM_QC}  ${PREBATCH_FLAG[*]:-}"
[[ -f "${SIM_QC}" ]] || { echo "ERROR: input not found: ${SIM_QC}"; exit 1; }
mkdir -p "${RUN}" "${SWEEP}/scores"

if [[ ! -f "${BANKSY_H5AD}" ]]; then
    echo "[1/3 banksy] building matrix (lambda=${LAM}, k_geom=${KG}, stagger=5, skip-umap, batch=${BATCH})"
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality "${MOD}" --input "${SIM_QC}" "${PREBATCH_FLAG[@]}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
else
    echo "[1/3 banksy] ${BANKSY_H5AD} exists, skipping build"
fi

if [[ ! -f "${RUN}/ari/ari_summary_${MOD}.json" ]]; then
    echo "[2/3 ari] resolution-matched Leiden ARI vs domain_true / cell_type_true"
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
        --modality "${MOD}" --packing-tag "${TAG}_${BATCH}" \
        --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"
else
    echo "[2/3 ari] ${RUN}/ari/ari_summary_${MOD}.json exists, skipping"
fi

echo "[3/3 composition] scoring (+ slice_id leakage ARI)"
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/misc/composition_recovery.py \
    --h5ad "${ARI_H5AD}" --tag "${MOD}_${BATCH}" \
    --out-dir "${SWEEP}/scores"

echo "[done] task ${i}  ->  ${RUN}"
