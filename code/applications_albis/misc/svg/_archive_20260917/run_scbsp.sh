#!/bin/bash
#SBATCH --job-name=scbsp_svg
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/svg/logs_scbsp/scbsp_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared      # uncomment / edit to match your cluster's partition name
#
# Figure 4A -- scBSP SVG identification, one array task per modality.
# CPU-only, no GPU needed (spot ~1s, bin16um/cell well under a minute on 4
# cores per scBSP's own benchmark scaling) -- run under sbatch anyway for
# logging/reproducibility consistency with the rest of the pipeline.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/svg/logs_scbsp

ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial"
PYTHON_BIN="${ENV_PREFIX}/bin/python"
export PYTHONUNBUFFERED=1

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

DATASETS=(bin16um spot cell)
DATASET="${DATASETS[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "Dataset: ${DATASET}"
echo "Command: ${PYTHON_BIN} sim_paper/code/applications_albis/svg/3D_scbsp.py --dataset ${DATASET} $*"
"${PYTHON_BIN}" sim_paper/code/applications_albis/svg/3D_scbsp.py --dataset "${DATASET}" "$@"

echo "Job finished: $(date)"
