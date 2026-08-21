#!/bin/bash
# Reproduces the bin (VisiumHD) noisy sweep: theta=200 (near-Poisson NB
# dispersion, up from the 2.0 default), to test whether raising theta alone
# can push the bin-level fit toward real breast_cancer_visium_hd (theta_hat
# ~1.23). Prior evidence (theta=2.0 -> theta_hat=0.05, theta=25 manuscript
# baseline -> theta_hat=0.07) suggests theta has little effect at the bin
# level -- pooling many spatially heterogeneous cells into one 8um bin is
# likely the dominant source of extra-Poisson variance there, not the
# per-molecule NB noise that theta controls. This run tests the ceiling.
#
# Submits the generate job, then the compare-vs-real job chained on it via
# --dependency=afterok so it runs automatically once generation finishes
# (~10 min for the bin modality).
#
# Usage: bash sim_paper/code/count_distribution/run_sweep_bin_theta200.sh
set -euo pipefail

cd /dcs04/hicks/data/Jan/sim_project

OUT_TAG=theta_200

GEN_JOB=$(sbatch --array=1-1 --parsable \
    sim_paper/code/data/run_generate_simulation_noisy.sh \
    --theta 200 --out-tag "${OUT_TAG}")
echo "Generate job (bin, theta=200): ${GEN_JOB}"

CMP_JOB=$(sbatch --parsable \
    --dependency=afterok:"${GEN_JOB}" \
    --export=ALL,OUT_TAG="${OUT_TAG}" \
    sim_paper/code/count_distribution/run_count_distribution_compare_noisy_visium_hd.sh)
echo "Compare job (vs breast_cancer_visium_hd): ${CMP_JOB}"
