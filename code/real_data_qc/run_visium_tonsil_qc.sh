#!/bin/bash
#SBATCH --job-name=visium_tonsil_qc
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/real_data_qc/logs/visium_tonsil_qc_%j.out
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# SpotSweeper QC for CytAssist FFPE Protein Expression Human Tonsil AddOns
# (probe-based, CytAssist; downloaded 2026-09-04 from
# cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Tonsil_AddOns/
# into data/real_data/tonsil_visium/ -- the processed spaceranger outputs
# (*_filtered_feature_bc_matrix.h5 + *_spatial.tar.gz), same bundle shape as
# lymph_node_visium; raw FASTQ/probe-set inputs not needed since these are
# already 10x's own spaceranger-processed matrix, 18,126 genes x 4,908 spots).
# Second 18k CytAssist probe reference for the Figure 2 "spot" comparison,
# alongside lymph_node_visium -- spot's real-reference set is being retuned
# jointly against both (breast_cancer_visium, 36k WTA, dropped per user
# decision 2026-09-04). Wraps visium_qc.py --sample tonsil -- see that file
# for what it does (SpotSweeper local-outlier filtering on
# total_counts/n_genes_by_counts/pct_counts_mt, drops flagged spots).
#
# Output under sim_paper/data/real_data_qc/tonsil_visium/:
#   tonsil_visium_qc.h5ad, qc_summary.json, qc_*.png
#
# Feeds into code/data/misc/sweep_spot_logmu_probe.sh and
# count_distribution/spot_vs_tonsil_visium.sh.
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
echo "Command: ${PYTHON_BIN} sim_paper/code/real_data_qc/visium_qc.py --sample tonsil $*"
"${PYTHON_BIN}" sim_paper/code/real_data_qc/visium_qc.py --sample tonsil "$@"

echo "Job finished: $(date)"
