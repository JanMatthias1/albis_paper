#!/bin/bash
#SBATCH --job-name=batch_slide_seeds
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/batch_slide_seeds_%A_%a.out
#SBATCH --array=0-71%8
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Seed-replication follow-up to cellbin_batch_sigma_slide.sh /
# batch_sigma_slide_fine.sh (batch_sigma_slide_domain_ari_final.png). Those
# swept batch_sigma at a SINGLE seed each; this adds 2 extra seeds (101, 202)
# per existing point across all three modalities, plus 2 new spot points
# beyond its canonical (0.4, 0.6), so the plot can show mean +/- spread
# instead of one line -- see gen_batch_sigma_slide_seed_tasks.py (rerun that
# to regenerate batch_sigma_slide_seed_tasks.tsv; this script just consumes
# it, one row per array task, same generate -> QC -> BANKSY+Harmony -> ARI ->
# composition_recovery pipeline as the two sweep scripts above).
#
# Requires FRESH simulations per task (the batch shift is baked into counts
# at generation time, same as the original sweeps).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/cellbin_batch_sigma_slide"
TASKS_TSV="sim_paper/code/clustering/cellbin_batch_sigma_slide/batch_sigma_slide_seed_tasks.tsv"

i="${SLURM_ARRAY_TASK_ID:-0}"
row=$((i + 2))  # +1 for 1-indexing, +1 for the header line
LINE="$(sed -n "${row}p" "${TASKS_TSV}")"
MOD="$(echo "${LINE}" | cut -f1)"
BS="$(echo "${LINE}" | cut -f2)"
SEED="$(echo "${LINE}" | cut -f3)"

FOLDER="${MOD}"
case "${MOD}" in
  cell)
    # 2026-09-18: synced to generate_strong_mix_cell.sh's authoritative config
    # (was stale at the OLD r=6000/log_mu=-2.3/k_geom=200 -- the exact combo
    # found broken/collapsed at this dataset's scale; see
    # project_figure3_banksy_domain_sweep memory 2026-09-18 for the full story).
    SPHERE_R_UM=2050
    LAM=0.5; KG=60
    GEN_FLAGS=(--n-cells 24207 --base-gene-lognormal -2.5 0.7 --theta 0.40 --theta-jitter 0.15 --strong-domain-mix)
    ;;
  bin)
    SPHERE_R_UM=2050
    LAM=0.5; KG=100
    GEN_FLAGS=(--bin-size-um 16 --base-gene-lognormal -2.5 0.7 --strong-domain-mix)
    # 2026-09-17: cellbin_batch_sigma_slide/ folder renamed "bin" -> "bin16um"
    # in the strong_mix consolidation -- --modality is still "bin" (only
    # cell/bin/spot are valid), but the ON-DISK folder is "bin16um".
    FOLDER="bin16um"
    ;;
  spot)
    # 2026-09-18: synced to generate_strong_mix_spot.sh's authoritative
    # dispersion (was stale at the bare generate_simulation_noisy.py defaults,
    # log_mu=-2.5/theta=2.0/jitter=1.0 -- predates Figure 2's spot retune to
    # log_mu=-2.0/theta=0.25/jitter=0.10, see [[reference_generate_noisy_theta_jitter_default]]).
    SPHERE_R_UM=2050
    LAM=0.1; KG=8
    GEN_FLAGS=(--base-gene-lognormal -2.0 0.7 --theta 0.25 --theta-jitter 0.10 --strong-domain-mix)
    ;;
esac

if [[ -n "${SEED}" ]]; then
    TAG="${MOD}_batch_slide_bs${BS}_seed${SEED}"
    RUN="${OUT_ROOT}/${FOLDER}/bs${BS}_seed${SEED}"
    SEED_FLAG=(--seed "${SEED}")
else
    TAG="${MOD}_batch_slide_bs${BS}"
    RUN="${OUT_ROOT}/${FOLDER}/bs${BS}"
    SEED_FLAG=()
fi
SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MOD}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MOD}_z_qc.h5ad"

echo "[task ${i}] modality=${MOD} batch_sigma=${BS} seed=${SEED:-<default>} lambda=${LAM} k_geom=${KG} -> ${RUN}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG}"
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MOD}" --sphere-r-um "${SPHERE_R_UM}" \
        "${GEN_FLAGS[@]}" \
        --batch-sigma "${BS}" \
        "${SEED_FLAG[@]}" \
        --out-tag "${TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${TAG}"
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/00_qc_filter.py \
        --modality "${MOD}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_${MOD}_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality "${MOD}" --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality "${MOD}" --packing-tag "cellbin_batch_slide_${TAG}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_${MOD}_z_ari_recovery.h5ad" --tag "${FOLDER}_$(basename "${RUN}")" \
    --out-dir "${OUT_ROOT}/scores"

# Render the manuscript panels from the completed embeddings and labels.
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/plot_banksy_results.py \
    --input "${BANKSY_H5AD}" \
    --cluster-input "${RUN}/ari/simulation_${MOD}_z_ari_recovery.h5ad"

echo "[done] task ${i}  ->  ${RUN}"
