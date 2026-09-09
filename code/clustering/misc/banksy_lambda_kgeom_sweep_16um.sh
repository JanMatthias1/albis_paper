#!/bin/bash
#SBATCH --job-name=banksy_lam_kg_16um
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/banksy_lam_kg_16um_%A_%a.out
#SBATCH --array=0-17%6
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# BANKSY lambda x k_geom sweep -- RE-RUN on the CURRENT shared strong-domain-mix
# tags (the byte-identical datasets Figure 4B STAGATE uses), replacing the
# 2026-08-25 configs the first sweep (banksy_lambda_kgeom_sweep.sh, job
# 35574059) used:
#     old bin : packing_pf0p04_bsigma05_strongdomainmix  (8um, log_mu 0.7, bsigma 0.5)
#     old cell: log_mu_-2.5_theta_0.25_strongdomainmix
# Here:
#     bin  -> 16um: packing_pf0p04_bin16um_log_mu_-2.5_bsigma07_strongmix
#     cell        : log_mu_-2.3_theta_0.40_jitter0.15_bsigma15_strongmix
# so Figure 3B (BANKSY domains) and Figure 4B (STAGATE domains) are on the
# same data and directly comparable. Both at the TUNED per-modality batch
# sigma (bin16 0.7, cell 1.5) -- BANKSY's PCA+Harmony does batch correction,
# so unlike STAGATE it may hold up; if the whole grid comes back ~0, retry on
# the ..._bsigma005_strongmix low-batch siblings.
#
# GRID (18 tasks): {cell, bin16um} x lambda {0.2, 0.5, 0.8} x k_geom {15, 60, 200}
#   -- the ORIGINAL sweep grid (lambda 0.2 cell-typing mode / 0.8 domain mode;
#   k_geom 15 the failed default, 60/200 push toward domain scale). 16um bins
#   are ~2x physically larger than the old 8um, so the same k_geom now spans
#   more tissue -- part of what this re-measures.
#   --stagger-scale 5 (no adjacent-slice leakage by construction),
#   --max-m 1, --skip-umap. composition_recovery.py reports ARI(clusters,
#   slice_id) as an explicit leakage flag: a run is a genuine domain hit only
#   if domain ARI is up AND slice_id-leakage ARI < ~0.3.
#
# Per task: (1) 01_build_banksy_matrix.py in sim-app-banksy env, (2)
# ari_vs_ground_truth.py + (3) composition_recovery.py --h5ad in
# albis-tutorial env (two-env handoff).
#
# Outputs under data/figure_3/banksy_lambda_kgeom_sweep_16um/:
#   <mod>/lam<L>_kg<K>/banksy_matrix/..._banksy_pca_harmony_qc.h5ad
#   <mod>/lam<L>_kg<K>/ari/ari_summary_<mod>.json
#   scores/composition_recovery_<mod>_lam<L>_kg<K>.json
# Tabulate with summary_banksy_lambda_kgeom.py (point SWEEP at this dir).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs

source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

SIM_DATA_DIR="sim_paper/data/figure_4/spatial_clustering/sim_data"
MODS=(cell bin)
LAMS=(0.2 0.5 0.8)
KGS=(15 60 200)

i="${SLURM_ARRAY_TASK_ID:-0}"
MOD="${MODS[$(( i / 9 ))]}"
rem=$(( i % 9 ))
LAM="${LAMS[$(( rem / 3 ))]}"
KG="${KGS[$(( rem % 3 ))]}"

case "${MOD}" in
  cell) SIM_TAG="log_mu_-2.3_theta_0.40_jitter0.15_bsigma15_strongmix"; MODFILE="cell" ;;
  bin)  SIM_TAG="packing_pf0p04_bin16um_log_mu_-2.5_bsigma07_strongmix"; MODFILE="bin"  ;;
esac
SIM_QC="${SIM_DATA_DIR}/${SIM_TAG}/simulation_${MODFILE}_z_qc.h5ad"
SWEEP="sim_paper/data/figure_3/banksy_lambda_kgeom_sweep_16um"
RUN="${SWEEP}/${MOD}/lam${LAM}_kg${KG}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_${MODFILE}_z_banksy_pca_harmony_qc.h5ad"
ARI_H5AD="${RUN}/ari/simulation_${MOD}_z_ari_recovery.h5ad"

echo "[task ${i}] MOD=${MOD} lambda=${LAM} k_geom=${KG}  tag=${SIM_TAG}"
echo "           input=${SIM_QC}"
[[ -f "${SIM_QC}" ]] || { echo "ERROR: input not found: ${SIM_QC}"; exit 1; }
mkdir -p "${RUN}" "${SWEEP}/scores"

if [[ ! -f "${BANKSY_H5AD}" ]]; then
    echo "[1/3 banksy] building matrix (lambda=${LAM}, k_geom=${KG}, stagger=5, skip-umap)"
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality "${MOD}" --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
else
    echo "[1/3 banksy] ${BANKSY_H5AD} exists, skipping build"
fi

if [[ ! -f "${RUN}/ari/ari_summary_${MOD}.json" ]]; then
    echo "[2/3 ari] resolution-matched Leiden ARI vs domain_true / cell_type_true"
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/ari_vs_ground_truth.py \
        --modality "${MOD}" --packing-tag "${SIM_TAG}" \
        --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"
else
    echo "[2/3 ari] ${RUN}/ari/ari_summary_${MOD}.json exists, skipping"
fi

echo "[3/3 composition] scoring (+ slice_id leakage ARI)"
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${ARI_H5AD}" --tag "${MOD}_lam${LAM}_kg${KG}" \
    --out-dir "${SWEEP}/scores"

echo "[done] task ${i}  ->  ${RUN}"
