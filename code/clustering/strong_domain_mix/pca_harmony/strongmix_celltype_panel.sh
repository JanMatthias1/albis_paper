#!/bin/bash
#SBATCH --job-name=strongmix_celltype_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/strongmix_celltype_panel_%A_%a.out
#SBATCH --array=0-17
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3, cell-type recovery panel on the STRONG domain-mix data -- the
# strong-mix counterpart of {cell,bin16um,spot}_celltype_panel.sh (which use
# the weak-mix Figure 2 data). Same pipeline: plain PCA->Harmony->Leiden,
# resolution-matched ARI vs cell_type_true, then UMAP + contingency plots.
#
# Inputs are the existing strong-mix QC files in
# data/figure_3/strong_domain_mix/batch_sigma_slide/<modality>/bs<X>/ (Figure 2 flags +
# --strong-domain-mix; nothing is regenerated here). Two batch levels per
# modality:
#   canonical   cell bs1.5 / bin16um bs0.7 / spot bs0.3  -- same batch_sigma
#               as the weak-mix panels, so weak vs strong mix differ only in
#               --strong-domain-mix
#   bs0         no batch effect -- same inputs as the strongmix rows of the
#               BANKSY lambda sweep (../banksy/)
# Outputs: data/figure_3/strong_domain_mix/pca_harmony/<modality>/bs<X>/
#   (spot's no-batch input folder is bs0.0; its output folder is bs0 like the others)
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

# task -> "<output name> <02/step03 --modality> <batch_sigma_slide dir>"
TASKS=(
    "cell cell cell/bs1.5"
    "cell cell cell/bs0"
    "bin16um bin bin16um/bs0.7"
    "bin16um bin bin16um/bs0"
    "spot spot spot/bs0.3"
    "spot spot spot/bs0.0"
    # 6-17: simulation seeds 101/202 (added 2026-09-25; submit --array=6-17 so
    # the seed-2025 runs above are not re-clustered)
    "cell cell cell/bs1.5_seed101"
    "cell cell cell/bs1.5_seed202"
    "cell cell cell/bs0_seed101"
    "cell cell cell/bs0_seed202"
    "bin16um bin bin16um/bs0.7_seed101"
    "bin16um bin bin16um/bs0.7_seed202"
    "bin16um bin bin16um/bs0_seed101"
    "bin16um bin bin16um/bs0_seed202"
    "spot spot spot/bs0.3_seed101"
    "spot spot spot/bs0.3_seed202"
    "spot spot spot/bs0_seed101"
    "spot spot spot/bs0_seed202"
)
read -r NAME MODALITY POINT <<< "${TASKS[${SLURM_ARRAY_TASK_ID}]}"

SLIDE_DIR="sim_paper/data/figure_3/strong_domain_mix/batch_sigma_slide/${POINT}"
QC_FILES=("${SLIDE_DIR}"/*_strongmix_bsigma*_qc.h5ad)
if [[ ${#QC_FILES[@]} -ne 1 || ! -f "${QC_FILES[0]}" ]]; then
    echo "[error] expected exactly one strong-mix QC h5ad in ${SLIDE_DIR}, found: ${QC_FILES[*]}" >&2
    exit 1
fi
SIM_QC="${QC_FILES[0]}"
OUT_POINT="$(basename "${POINT}")"
[[ "${OUT_POINT}" == bs0.0 ]] && OUT_POINT=bs0
CLUSTER_ROOT="sim_paper/data/figure_3/strong_domain_mix/pca_harmony/${NAME}/${OUT_POINT}"
PCA_H5AD="${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad"
echo "[input] ${SIM_QC} -> ${CLUSTER_ROOT}"

if [[ ! -f "${PCA_H5AD}" ]]; then
    echo "[pca_harmony] running"
    "${PYTHON_BIN}" sim_paper/code/clustering/step01_pca_harmony.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --no-umap-sample \
        --output "${PCA_H5AD}"
fi

echo "[ari] resolution-matched ARI recovery"
"${PYTHON_BIN}" sim_paper/code/clustering/step02_leiden_resolution_sweep.py \
    --modality "${MODALITY}" \
    --input "${PCA_H5AD}" \
    --output-dir "${CLUSTER_ROOT}/ari_recovery_qc"

# Same cell_type_true-matched resolution as the weak-mix panels (see
# cell_celltype_panel.sh for why).
RESOLUTION=$("${PYTHON_BIN}" -c "
import json
with open('${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json') as f:
    summary = json.load(f)
print(next(r['resolution'] for r in summary if r['ground_truth'] == 'cell_type_true'))
")
echo "[leiden] cell_type_true-matched-resolution qualitative plots (resolution=${RESOLUTION})"
"${PYTHON_BIN}" sim_paper/code/clustering/step03_cluster_and_plot.py \
    --modality "${MODALITY}" \
    --input "${PCA_H5AD}" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_qc_celltype_matched" \
    --pipeline pca_harmony --resolution "${RESOLUTION}"

echo "[done] cell_type_true ARI -> ${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json"
