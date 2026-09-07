#!/bin/bash
#SBATCH --job-name=regen_qc_postbatch
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/regen_qc_postbatch_%j.out
#SBATCH --time=02:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# 2026-09-01: regenerate the 4 canonical Figure 2 QC h5ads after 00_qc_filter.py
# was switched from filtering on counts_pre_batch to the post-batch matrix
# (adata.X). The old _qc.h5ad let near-empty bins pass the >=3-genes floor on
# their pre-batch counts, then the batch noise zeroed them -> detached-island
# artifact in the Fig 3 UMAPs and a near-empty tail in the Fig 2 comparisons.
# Overwrites data/figure_2/<tag>/simulation_<mod>_z_qc.h5ad in place.

set -euo pipefail
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

run_qc () {
    local tag="$1" mod="$2"
    local raw="sim_paper/data/figure_2/${tag}/simulation_${mod}_z.h5ad"
    local qc="sim_paper/data/figure_2/${tag}/simulation_${mod}_z_qc.h5ad"
    echo "=== ${tag} (${mod}) ==="
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
        --modality "${mod}" --input "${raw}" --output "${qc}"
}

run_qc "packing_pf0p04_log_mu_0.0_bsigma08"            bin
run_qc "packing_pf0p04_bin16um_log_mu_-2.5_bsigma07"   bin
run_qc "packing_pf0p04_log_mu_-2.5_bsigma03"           spot
run_qc "log_mu_-2.3_theta_0.40_jitter0.15_bsigma15"    cell

echo "[done] all 4 QC h5ads regenerated on post-batch counts"
