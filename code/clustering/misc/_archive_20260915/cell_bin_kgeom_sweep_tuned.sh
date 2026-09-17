#!/bin/bash
#SBATCH --job-name=cellbin_kgeom_sweep
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/cellbin_kgeom_sweep_%A_%a.out
#SBATCH --array=0-15%8
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Follow-up to spot_kgeom_sweep_tuned.sh's finding: at the real tuned batch
# effect, spot's domain recovery collapses at k_geom>=15 into slice leakage,
# but k_geom=6 (lambda=0.2) gets a clean, if modest, ARI=0.241 through.
# Cell/bin's original sweep (job 35589174) never tried below k_geom=15 --
# this tests whether the same small-k_geom effect rescues them too, or
# whether spot's win was specific to its much sparser density.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

SIM_DATA_DIR="sim_paper/data/figure_4/spatial_clustering/sim_data"
OUT_ROOT="sim_paper/data/figure_3/cellbin_kgeom_sweep_tuned"

# 8 (lambda, k_geom) combos x 2 modalities = 16 tasks
LAMBDAS=(0.2 0.2 0.2 0.2 0.5 0.5 0.5 0.5)
KGEOMS=(4 6 10 15 4 6 10 15)
MODS=(cell bin)
TAGS=(
  "log_mu_-2.3_theta_0.40_jitter0.15_bsigma15_strongmix"
  "packing_pf0p04_bin16um_log_mu_-2.5_bsigma07_strongmix"
)

i="${SLURM_ARRAY_TASK_ID:-0}"
COMBO_IDX=$(( i % 8 ))
MOD_IDX=$(( i / 8 ))
LAM="${LAMBDAS[$COMBO_IDX]}"
KG="${KGEOMS[$COMBO_IDX]}"
MOD="${MODS[$MOD_IDX]}"
TAG_DATA="${TAGS[$MOD_IDX]}"

SIM_QC="${SIM_DATA_DIR}/${TAG_DATA}/simulation_${MOD}_z_qc.h5ad"
TAG="lambda${LAM}_kgeom${KG}"
RUN="${OUT_ROOT}/${MOD}/${TAG}"

echo "[task ${i}] modality=${MOD} lambda=${LAM} k_geom=${KG} (tuned/real batch)"
mkdir -p "${RUN}"

BANKSY_H5AD="${RUN}/banksy_matrix/simulation_${MOD}_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality "${MOD}" --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/ari_vs_ground_truth.py \
    --modality "${MOD}" --packing-tag "cellbin_kgeom_sweep_${MOD}_${TAG}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_${MOD}_z_ari_recovery.h5ad" --tag "${MOD}_${TAG}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] task ${i}  ->  ${RUN}"
