#!/bin/bash
#SBATCH --job-name=gen_cell_matched_scale
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/gen_cell_matched_scale_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# One-off: regenerate the "cell" modality at the SAME sphere_r_um=2050 disc
# bin16um/spot already use (packing_pf0p04 family), instead of cell's
# manuscript-default sphere_r_um=6000. All other params are held identical
# to the existing canonical cell tag (log_mu_-2.3_theta_0.40_jitter0.15_
# bsigma15, uns['sim_params'] confirmed) -- only --sphere-r-um changes.
#
# Why: cross_tech_stair.py (Figure 4E) currently rescales cell's coords onto
# the bin/spot disc post-hoc (factor ~0.343) before STAIR alignment, because
# a nearest-neighbour check showed all three modalities already carry the
# SAME domain_true field, just at different absolute scales -- the domain
# layout appears to be defined relative to sphere_r_um, not by an absolute
# um coordinate. Generating cell natively at sphere_r_um=2050 should put it
# in the same physical frame as bin16um/spot from the start, removing that
# rescale step and the packing-fraction confound (native cell packs 600k
# cells at ~0.16%; this run packs the same 600k into the bin/spot disc, so
# ~4% -- like_packing_pf0p04 -- much closer to the "real tissue is nearly
# fully packed" note in generate_simulation_noisy.py's header comment).
#
# NOT written to data/figure_2/ -- this is exploratory for Figure 4E only.
# Stays under data/noisy/<tag>/ until the cross-tech alignment result on it
# is checked; promote manually (mv to data/figure_2/) only if adopted.
#
# Usage: sbatch generate_cell_matched_scale.sh

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

TAG="packing_pf0p04_cell_log_mu_-2.3_theta_0.40_jitter0.15_bsigma15"
NOISY_DIR="sim_paper/data/noisy/${TAG}"

echo "[config] modality=cell log_mu=-2.3 theta=0.40 jitter=0.15 batch_sigma=1.5 sphere_r_um=2050 tag=${TAG}"

if [[ -f "${NOISY_DIR}/simulation_cell_z_qc.h5ad" ]]; then
    echo "[skip] ${NOISY_DIR} already complete"
    exit 0
fi

GEN_ARGS=(--modality cell --base-gene-lognormal -2.3 0.7 --batch-sigma 1.5 \
          --theta 0.40 --theta-jitter 0.15 --sphere-r-um 2050 --out-tag "${TAG}")

if [[ ! -f "${NOISY_DIR}/simulation_cell_z.h5ad" ]]; then
    echo "[generate] ${TAG}"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py "${GEN_ARGS[@]}"
fi

if [[ ! -f "${NOISY_DIR}/simulation_cell_z_qc.h5ad" ]]; then
    echo "[qc] ${TAG}"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality cell --packing-tag "${TAG}"
fi

echo "[done] ${TAG} -> ${NOISY_DIR}"
