#!/bin/bash
#SBATCH --job-name=visium_lymph_qc
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/real_data_qc/logs/visium_lymph_qc_%j.out
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# SpotSweeper QC for the Visium V2 Human Lymph Node sample (probe-based
# CytAssist FFPE; downloaded from
# cf.10xgenomics.com/samples/spatial-exp/3.1.3/Visium_V2_Human_Lymph_Node/
# into data/real_data/lymph_node_visium/) -- one of the two real Visium
# references for the Figure 2 "spot" comparison, alongside tonsil_visium. Wraps
# visium_qc.py --sample lymph_node -- see that file for what it does
# (SpotSweeper local-outlier filtering on total_counts/n_genes_by_counts/
# pct_counts_mt, drops flagged spots).
#
# Output under sim_paper/data/real_data_qc/lymph_node_visium/:
#   lymph_node_visium_qc.h5ad, qc_summary.json, qc_*.png
#
# Feeds into count_distribution/smaller_sphere/spot_vs_lymph_node_visium.sh.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/real_data_qc/logs

# visium_qc.py only needs scanpy/spotsweeper/matplotlib -- same env
# count_distribution.py uses. Absolute path (not `command -v python` after
# conda activate) -- on a busy node that can silently resolve to a different
# env's interpreter still ahead on PATH (see count_distribution/_env.sh).
ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial"
PYTHON_BIN="${ENV_PREFIX}/bin/python"

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"
echo "Command: ${PYTHON_BIN} sim_paper/code/real_data_qc/visium_qc.py --sample lymph_node $*"
"${PYTHON_BIN}" sim_paper/code/real_data_qc/visium_qc.py --sample lymph_node "$@"

echo "Job finished: $(date)"
