#!/bin/bash
#SBATCH --job-name=strong_mix_cell
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_mix/logs/strong_mix_cell_%A_%a.out
#SBATCH --array=0-8%4
#SBATCH --time=06:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Cell's strong-domain-mix batch_sigma slide, consolidated 2026-09-17 from
# the scattered points built across code/clustering/misc/
# cellbin_batch_sigma_slide.sh (bs0.5/1.0) + batch_sigma_slide_fine.sh
# (bs0.1/0.2/0.3/0.4) into one self-contained array script, output moved
# from data/figure_3/cellbin_batch_sigma_slide/cell/ to data/strong_mix/cell/.
# Config: BANKSY lambda=0.5 (cell's own Phase-1 ceiling, unchanged), k_geom=60
# -- RE-TUNED 2026-09-18, see below (was 200, cell's old Phase-1 ceiling from
# a ~25x-larger cell count -- broke completely at this dataset's new scale).
#
# 2026-09-17 (same day, later): sphere_r_um switched from cell's old native
# 6000 to 2050 + n_cells added (24207) -- matching bin16um/spot's
# smaller_sphere disc, per user request ("generate cell strong mix with the
# same diameter as well"). Before this, spot/bin16um shared the 2050 disc
# while cell alone sat on its own 6000 sphere with an implicit default
# n_cells=600000 (generate_simulation_noisy.py does NOT auto-derive n_cells
# from sphere_r_um -- would have silently stayed 600000 even if only the
# radius were changed). 24207 is the exact n_cells Figure 2's smaller_sphere
# cell tag already uses (chosen from a dedicated fidelity sweep across
# 24207/75000/200000/600000 at r=2050 -- 24207 scored best; see
# project_figure2_smaller_larger_sphere memory), reused here rather than
# re-deriving a new value. Old r=6000/n=600000 results (2026-09-17, same-day)
# archived to cellbin_batch_sigma_slide/cell_pre_diameter_match_20260917/ and
# scores/_pre_diameter_match_20260917/. Resources dropped 320G/18h -> 96G/6h
# to match the ~25x smaller per-point cell count (was sized for 600k cells).
#
# 2026-09-17 (same day, later still): dispersion log_mu -2.3 -> -2.5, matching
# Figure 2's own retune -- an actual joint log_mu x theta grid search against
# both Xenium refs (code/data/misc/sweep_cell_logmu_theta_joint_xenium.sh, job
# 35759008), since cell had never gotten spot's joint-sweep treatment before
# (see project_figure4c_alignment_4modality memory for the full sweep table
# and Figure 2's before/after comparison -- total_counts fidelity ratio
# 1.38->1.04, theta_hat/genes_per_cell shift slightly further from real as
# the accepted tradeoff). theta/jitter unchanged (0.40/0.15, already
# artifact-free).
#
# Grid EXPANDED from 6 to 9 points (0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0,
# 1.5): bs=0 and the canonical bs=1.5 used to live separately under
# data/figure_3/banksy_batch_compare/cell/{prebatch,tuned}/ at the OLD
# r=6000/n=600000/log_mu=-2.3 config (confirmed directly from their h5ad's
# obsm['spatial'], max radius ~5989um) -- diameter- AND dispersion-
# inconsistent with the rest of this curve. bs=0.05 was a newer addition
# (low-batch reference for Figure 4B's STAGATE) at the same stale config.
# All 3 now generated here instead, at the matched geometry+dispersion, so
# this script is the single source of truth for cell's entire strong-mix
# curve. IMPORTANT: bs=1.5 and bs=0.05 are read directly by
# code/applications_albis/spatial_clustering/3D_stagate.py (Figure 4B,
# "cell_strongmix"/"cell_strongmix_lowbatch" datasets) -- regenerating them
# here changes Figure 4B's input too; that panel needs a rerun after this
# lands (see project_figure4c_alignment_4modality memory). Old bs0/bs0.05/
# bs1.5 archived alongside the other 6 points (see above).
#
# 2026-09-18: k_geom RE-TUNED 200 -> 60 after running the above 9-point grid
# and finding bs=0.2 through bs=1.5 (6 of 9 points) came back BYTE-IDENTICAL
# (ari=0.0017370357583088834 to 16 decimal places, slice_id_leakage_ari=1.0
# for every one) -- impossible by chance from 6 independently-generated
# datasets, and confirmed as a k_geom-driven collapse: k_geom=200 was tuned
# for the OLD ~60,000-cells/slice scale; at the new ~2,400/slice scale (as
# few as 683 in the thinnest slice) it smooths over 7-29% of an entire
# slice per cell. Even bs=0 (no batch effect) only reached domain ARI 0.077
# at k_geom=200, vs 0.784 in the old, correctly-scaled dataset -- this
# wasn't a batch-robustness finding, k_geom=200 simply doesn't work at this
# scale. Dedicated re-tune (code/clustering/misc/sweep_cell_kgeom_r2050.sh,
# job 35786358, k_geom in {15,30,60,100,150,200} at bs=0, lambda=0.5 held):
#   k_geom=60   domain_ARI=0.5404  leakage=0.0524   <- winner
#   k_geom=30   domain_ARI=0.5072  leakage=0.0090
#   k_geom=100  domain_ARI=0.3550  leakage=0.0514
#   k_geom=15   domain_ARI=0.3447  leakage=0.0023
#   k_geom=150  domain_ARI=0.1242  leakage=0.4421   <- leakage blow-up starts
#   k_geom=200  domain_ARI=0.1062  leakage=0.5333   <- confirms the collapse
# 60 matches a naive proportional-scaling estimate (200 * 2050/6000 = 68)
# reasonably well -- consistent with the mechanism being that
# STRONG_DOMAIN_MIX domains scale in physical size with sphere_r_um, so a
# fixed k_geom neighborhood now spans a much larger fraction of a domain's
# diameter than before. Old k_geom=200 results (2026-09-18) archived to
# cellbin_batch_sigma_slide/cell_pre_kgeom_retune_20260918/ and
# scores/_pre_kgeom_retune_20260918/. Only banksy_matrix/ari/composition_recovery
# need rerunning at the new k_geom -- raw sim data (data/noisy/cell_strong_mix_bs*/)
# is k_geom-independent and untouched, so this rerun is cheap.
#
# Full pipeline per point: generate (strong-domain-mix) -> QC ->
# BANKSY+Harmony -> 02_leiden_resolution_sweep.py -> composition_recovery.py
# (ARI + slice_id leakage).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_mix/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/cellbin_batch_sigma_slide"
LAM=0.5
KG=60

BATCH_SIGMAS=(0 0.05 0.1 0.2 0.3 0.4 0.5 1.0 1.5)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="cell_strong_mix_bs${BS}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_cell_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_cell_z_qc.h5ad"
RUN="${OUT_ROOT}/cell/bs${BS}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] cell batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality cell --sphere-r-um 2050 --n-cells 24207 \
        --base-gene-lognormal -2.5 0.7 --theta 0.40 --theta-jitter 0.15 \
        --strong-domain-mix --batch-sigma "${BS}" \
        --out-tag "${TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/00_qc_filter.py \
        --modality cell --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_cell_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality cell --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality cell --packing-tag "strong_mix_cell_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_cell_z_ari_recovery.h5ad" --tag "cell_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] task ${SLURM_ARRAY_TASK_ID:-0} -> ${RUN}"
