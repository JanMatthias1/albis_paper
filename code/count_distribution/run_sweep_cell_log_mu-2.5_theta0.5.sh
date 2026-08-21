#!/bin/bash
# Reproduces the cell-modality noisy sweep: base_gene_lognormal log_mu=-2.5
# (same as the log_mu_-2.5 run) combined with theta=0.5 (down from the 2.0
# default), to test whether lowering theta further closes the remaining
# overdispersion gap vs real Xenium lung data (sim theta_hat ~0.76 vs real
# ~0.15-0.16 at theta=2.0).
#
# Submits the generate job, then the compare-vs-real jobs (both lung slices)
# chained on it via --dependency=afterok so they run automatically once
# generation finishes (~20 min for 600k cells).
#
# Usage: bash sim_paper/code/count_distribution/run_sweep_cell_log_mu-2.5_theta0.5.sh
set -euo pipefail

cd /dcs04/hicks/data/Jan/sim_project

OUT_TAG=log_mu_-2.5_theta_0.5

GEN_JOB=$(sbatch --array=2-2 --parsable \
    sim_paper/code/data/run_generate_simulation_noisy.sh \
    --base-gene-lognormal -2.5 0.7 --theta 0.5 --out-tag "${OUT_TAG}")
echo "Generate job (cell, theta=0.5): ${GEN_JOB}"

CMP_JOB=$(sbatch --array=0-1 --parsable \
    --dependency=afterok:"${GEN_JOB}" \
    --export=ALL,OUT_TAG="${OUT_TAG}" \
    sim_paper/code/count_distribution/run_count_distribution_compare_noisy.sh)
echo "Compare job (vs non_diseased_lung, lung_cancer): ${CMP_JOB}"
