#!/bin/bash
#SBATCH --job-name=sweep_cell_kgeom_r2050
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/sweep_cell_kgeom_r2050_%A_%a.out
#SBATCH --array=0-5
#SBATCH --time=04:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# 2026-09-18: cell's BANKSY k_geom re-tune, forced by the same-day diameter
# fix (strong_mix/generate_strong_mix_cell.sh, r=6000/n=600000 ->
# r=2050/n=24207, matching bin16um/spot's smaller_sphere disc). k_geom=200
# was cell's "Phase-1 ceiling" from the OLD ~60,000-cells/slice scale; at the
# new ~2,400/slice scale (as few as 683 in the thinnest slice) it collapsed
# domain recovery completely -- confirmed by running the full 9-point
# strong-mix grid at k_geom=200: bs=0.2 through bs=1.5 (6 of 9 points) came
# back with BYTE-IDENTICAL ARI (0.0017370357583088834 to 16 decimal places)
# and slice_id_leakage_ari=1.0 for every one of them. That is not batch-effect
# behavior -- since batch_sigma changes only expression noise, not the
# ground-truth domain_true/cell_type_true assignment, an ARI that is
# constant across 6 different batch_sigma values means clustering converged
# to the SAME degenerate partition (== slice_id) regardless of batch, which
# only k_geom-driven over-smoothing explains. Even bs=0 (no batch effect at
# all) only reached domain ARI 0.077, vs 0.784 in the old, correctly-scaled
# dataset -- so this isn't a batch-robustness question, k_geom=200 is simply
# too large a neighborhood at this scale (200 neighbors is 7-29% of a
# slice's total cells, vs ~0.33% at the old scale).
#
# Mechanism (most likely): STRONG_DOMAIN_MIX domains scale in physical size
# with sphere_r_um (smaller sphere -> smaller domains in absolute microns);
# k_geom's effective smoothing radius depends on local cell density, which
# IS held constant by design (n_cells scaled by (r_new/r_old)**3 specifically
# to preserve packing fraction -- see project_figure2_smaller_larger_sphere
# memory) -- so a fixed k_geom now spans a much larger FRACTION of a
# domain's diameter than before, smoothing across domain boundaries it
# previously stayed within. Naive proportional scaling (200 * 2050/6000 =
# 68) motivates this grid's range.
#
# Grid: k_geom in {15, 30, 60, 100, 150, 200} (lambda=0.5 held fixed --
# cell's existing choice, not implicated by the scale change; this sweep
# only re-checks the k_geom axis). Run against the bs=0 (no-batch) point
# from the just-generated 9-point strong-mix grid -- isolates "can BANKSY
# find real domain structure at this scale at all" from the separate
# question of batch-robustness (matches bin16um_lambda_kgeom_sweep_16um.sh's
# own fallback: "if the whole grid comes back ~0, retry on ..._lowbatch").
# Reuses the ALREADY-GENERATED QC'd h5ad (data/noisy/cell_strong_mix_bs0/) --
# only the BANKSY+clustering steps are swept, not resimulated.
#
# Outputs under data/figure_3/cellbin_batch_sigma_slide/misc/cell_kgeom_sweep_r2050/:
#   kg<K>/banksy_matrix/..._banksy_pca_harmony_qc.h5ad
#   kg<K>/ari/ari_summary_cell.json
#   scores/composition_recovery_cell_kg<K>.json
# Tabulate with summary_cell_kgeom_r2050.sh (companion script).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

SIM_QC="sim_paper/data/noisy/cell_strong_mix_bs0/simulation_cell_z_qc.h5ad"
SWEEP="sim_paper/data/figure_3/cellbin_batch_sigma_slide/misc/cell_kgeom_sweep_r2050"
LAM=0.5

KGS=(15 30 60 100 150 200)
KG="${KGS[${SLURM_ARRAY_TASK_ID:-0}]}"
RUN="${SWEEP}/kg${KG}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_cell_z_banksy_pca_harmony_qc.h5ad"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] cell k_geom=${KG} lambda=${LAM} (bs=0, r=2050/n=24207)"
[[ -f "${SIM_QC}" ]] || { echo "ERROR: input not found: ${SIM_QC}"; exit 1; }
mkdir -p "${RUN}" "${SWEEP}/scores"

if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality cell --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality cell --packing-tag "cell_kgeom_sweep_r2050_kg${KG}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_cell_z_ari_recovery.h5ad" --tag "cell_kg${KG}" \
    --out-dir "${SWEEP}/scores"

echo "[done] task ${SLURM_ARRAY_TASK_ID:-0} -> ${RUN}"
