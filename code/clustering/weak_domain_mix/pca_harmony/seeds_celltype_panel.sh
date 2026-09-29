#!/bin/bash
#SBATCH --job-name=weakmix_seeds_celltype_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/weak_domain_mix/logs/seeds_celltype_panel_%A_%a.out
#SBATCH --array=0-11
#SBATCH --time=06:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Plain PCA->Harmony->Leiden on the weak domain-mix simulation seeds 101/202
# (../generate/generate_weakmix_seeds.sh), at the usual batch_sigma (Figure 2
# config) and at batch_sigma 0: resolution-matched ARI vs cell_type_true and
# domain_true, then UMAP + contingency plots. Seed 2025 is covered by
# {cell,bin16um,spot}_celltype_panel.sh (usual sigma) and bs0_celltype_panel.sh.
# Gives every Figure 3 bar three seeds (mean ± SD).
#
# Task = modality (3) x batch (usual, 0) x seed (101, 202), same order as the
# generator. Output: data/figure_3/weak_domain_mix/pca_harmony/<modality>/bs<sigma>_seed<n>/
# DRY_RUN=1 prints the commands instead of running them.
set -euo pipefail

source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project
RUN=(); [[ -n "${DRY_RUN:-}" ]] && RUN=(echo)

TASK="${SLURM_ARRAY_TASK_ID:?run as an array task}"
MODS=(cell bin16um spot); BATCHES=(usual 0); SEEDS=(101 202)
NAME="${MODS[$((TASK / 4))]}"; BATCH="${BATCHES[$(((TASK / 2) % 2))]}"; SEED="${SEEDS[$((TASK % 2))]}"

UMAP_FLAGS=(--no-umap-sample)   # full-data UMAP; bin16um overrides below
case "${NAME}" in
    cell) MOD=cell; SIGMA=1.5; TAG="log_mu_-2.5_theta_0.40_jitter0.15_bsigma15" ;;
    bin16um) MOD=bin; SIGMA=0.7; TAG="packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma07"
          # full UMAP on ~366k bins takes >5h; subsample it (plots only --
          # PCA/Harmony/Leiden/ARI use every bin)
          UMAP_FLAGS=() ;;
    spot) MOD=spot; SIGMA=0.3; TAG="packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_bsigma03" ;;
esac
if [[ "${BATCH}" == usual ]]; then
    SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/${TAG}_seed${SEED}/simulation_${MOD}_z_qc.h5ad"
    POINT="bs${SIGMA}"
else
    QC_FILES=(sim_paper/data/figure_3/weak_domain_mix/data/${NAME}_seed${SEED}/*_bsigma0_seed${SEED}_qc.h5ad)
    SIM_QC="${QC_FILES[0]}"
    POINT="bs0"
fi
CLUSTER_ROOT="sim_paper/data/figure_3/weak_domain_mix/pca_harmony/${NAME}/${POINT}_seed${SEED}"
PCA_H5AD="${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad"
ARI_JSON="${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MOD}.json"
echo "[task ${TASK}] ${NAME} batch=${BATCH} seed=${SEED}: ${SIM_QC} -> ${CLUSTER_ROOT}"
[[ -n "${DRY_RUN:-}" || -f "${SIM_QC}" ]] || { echo "[error] input missing: ${SIM_QC}" >&2; exit 1; }

if [[ ! -f "${PCA_H5AD}" ]]; then
    "${RUN[@]}" "${PYTHON_BIN}" sim_paper/code/clustering/step01_pca_harmony.py \
        --modality "${MOD}" --input "${SIM_QC}" "${UMAP_FLAGS[@]}" --output "${PCA_H5AD}"
fi
if [[ ! -f "${ARI_JSON}" ]]; then
    "${RUN[@]}" "${PYTHON_BIN}" sim_paper/code/clustering/step02_leiden_resolution_sweep.py \
        --modality "${MOD}" --input "${PCA_H5AD}" --output-dir "${CLUSTER_ROOT}/ari_recovery_qc"
fi
[[ -n "${DRY_RUN:-}" ]] && { echo "[dry-run] then step03_cluster_and_plot.py at the cell_type_true-matched resolution"; exit 0; }

# Same cell_type_true-matched resolution as the other panels (see
# cell_celltype_panel.sh for why).
RESOLUTION=$("${PYTHON_BIN}" -c "
import json
with open('${ARI_JSON}') as f:
    summary = json.load(f)
print(next(r['resolution'] for r in summary if r['ground_truth'] == 'cell_type_true'))
")
"${PYTHON_BIN}" sim_paper/code/clustering/step03_cluster_and_plot.py \
    --modality "${MOD}" \
    --input "${PCA_H5AD}" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_qc_celltype_matched" \
    --pipeline pca_harmony --resolution "${RESOLUTION}"
echo "[done] ${ARI_JSON}"
