#!/bin/bash
#SBATCH --job-name=banksy_lam_kg_sweep
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/banksy_lam_kg_sweep_%A_%a.out
#SBATCH --array=0-17%6
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# BANKSY lambda x k_geom sweep on the STRONG-domain-mix bin + cell datasets.
#
# WHY: on strong-mix data BANKSY (lambda=0.8, k_geom=15, stagger=1.5) helps spot
# (domain ARI 0.37 -> 0.61) but for bin/cell it collapses BOTH domain and
# cell-type recovery to ~0, and the k_geom=30-100 variants that looked better
# were a slice-leakage artifact (staggered coords -> neighbor search crosses
# into the adjacent slice -> clusters track slice_id, not domain). Question
# (PI, 2026-09-06): is bin/cell just a BANKSY parameter-tuning gap?
#
# GRID (18 tasks): {cell, bin} x lambda {0.2, 0.5, 0.8} x k_geom {15, 60, 200}.
#   lambda 0.2 = cell-typing mode (mostly own expression), 0.8 = domain mode.
#   k_geom  15 = the failed default; 60/200 push the spatial neighbourhood
#           toward domain scale (still small vs a ~1-2 mm domain at bin/cell
#           density -- part of what this measures).
#   --stagger-scale 5 (up from the default 1.5) so the larger k_geom neighbour
#   searches cannot reach the adjacent staggered slice -> no leakage by
#   construction; composition_recovery.py also reports ARI(clusters, slice_id)
#   as an explicit leakage flag.
#   --skip-umap: only X_pca_harmony is scored; UMAP on 600-716k obs would
#   dominate a sweep this size.
#
# Per task: (1) 01_build_banksy_matrix.py in the sim-app-banksy env, (2)
# ari_vs_ground_truth.py + (3) composition_recovery.py --h5ad in the
# albis-tutorial env (two-env handoff, same pattern as *_domain_panel.sh).
#
# Outputs under data/figure_3/banksy_lambda_kgeom_sweep/:
#   <mod>/lam<L>_kg<K>/banksy_matrix/simulation_<mod>_z_banksy_pca_harmony_qc.h5ad
#   <mod>/lam<L>_kg<K>/ari/ari_summary_<mod>.json
#   scores/composition_recovery_<mod>_lam<L>_kg<K>.json
# Tabulate with: code/clustering/misc/summary_banksy_lambda_kgeom.sh
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs

source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

MODS=(cell bin)
LAMS=(0.2 0.5 0.8)
KGS=(15 60 200)

i="${SLURM_ARRAY_TASK_ID:-0}"
MOD="${MODS[$(( i / 9 ))]}"
rem=$(( i % 9 ))
LAM="${LAMS[$(( rem / 3 ))]}"
KG="${KGS[$(( rem % 3 ))]}"

case "${MOD}" in
  cell) SIM_TAG="log_mu_-2.5_theta_0.25_strongdomainmix"; MODFILE="cell" ;;
  bin)  SIM_TAG="packing_pf0p04_bsigma05_strongdomainmix"; MODFILE="bin" ;;
esac
SIM_QC="sim_paper/data/noisy/${SIM_TAG}/simulation_${MODFILE}_z_qc.h5ad"
SWEEP="sim_paper/data/figure_3/banksy_lambda_kgeom_sweep"
RUN="${SWEEP}/${MOD}/lam${LAM}_kg${KG}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_${MODFILE}_z_banksy_pca_harmony_qc.h5ad"
ARI_H5AD="${RUN}/ari/simulation_${MOD}_z_ari_recovery.h5ad"

echo "[task ${i}] MOD=${MOD} lambda=${LAM} k_geom=${KG}  tag=${SIM_TAG}"
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
