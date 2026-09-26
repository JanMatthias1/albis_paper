#!/usr/bin/env bash
# One-shot follow-up to job 35789683 (generate_strong_mix_bin16um.sh, array
# task 0 = the bin16um bs0 "poisson baseline" regen). Submit this WITH
# --dependency=afterok:35789683_0 so it only runs once that job has landed
# fresh ari/banksy_matrix output for data/figure_3/strong_domain_mix/batch_sigma_slide/
# bin16um/bs0/ -- both scripts below read that folder (plot_domain_vs_celltype.py
# directly, plot_batch_sigma_slide_final.py via its auto-discovery over all
# bs<value> folders) and were left pointed at "bs0" (not the pre-regen
# archive) in anticipation of this landing. See
# project_figure3_banksy_domain_sweep.md's 2026-09-18 entry for the full
# thread.
#SBATCH --job-name=replot_bin16um_bs0
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/ari_recovery_summary/logs/replot_bin16um_bs0_%j.out
#SBATCH --time=00:20:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=shared
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/ari_recovery_summary/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

echo "[replot] plot_domain_vs_celltype.py"
"${PYTHON_BIN}" sim_paper/code/clustering/ari_recovery_summary/plot_domain_vs_celltype.py

echo "[replot] plot_batch_sigma_slide_final.py"
"${PYTHON_BIN}" sim_paper/code/clustering/strong_domain_mix/batch_sigma_slide/plot_batch_sigma_slide_final.py

echo "[done] both figures regenerated against bin16um/bs0's new poisson-baseline data"
