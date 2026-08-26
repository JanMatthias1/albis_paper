#!/bin/bash
#SBATCH --job-name=visium_hd_qc_16um
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/real_data_qc/logs/visium_hd_qc_16um_%A_%a.out
#SBATCH --array=0-1
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Real Visium HD 16um-bin references, same SpotSweeper QC pipeline already
# used for the 8um references (visium_hd_qc.py --bin-size-um defaults to 8;
# the raw 10x binned_outputs/ already ships square_016um alongside
# square_008um, so this needs no new download). For the bin_size_um=16
# sweep comparing against 8um's established config (packing_pf0p04_log_mu_0.0),
# see figure.md Figure 3 bin BANKSY domain-recovery discussion.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/real_data_qc/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/real_data_qc/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SAMPLES=(breast_cancer human_pancreas)
SAMPLE="${SAMPLES[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "Sample: ${SAMPLE}"
"${PYTHON_BIN}" sim_paper/code/real_data_qc/visium_hd_qc.py \
    --sample "${SAMPLE}" --bin-size-um 16 \
    --output-dir "sim_paper/data/real_data_qc/${SAMPLE}_visium_hd_16um"
