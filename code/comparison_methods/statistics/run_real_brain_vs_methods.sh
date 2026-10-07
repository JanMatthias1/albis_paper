#!/usr/bin/env bash
#SBATCH --job-name=real_brain_vs_methods
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/comparison_methods/statistics/logs/real_brain_vs_methods_%j.out
# Each QC'd Xenium brain sample (code/real_data_qc/run_brain_xenium_qc.sh) vs the Figure 5B cell
# inputs (ALBIS, SPIDER, scCube; slice 4, QC'd), simulated panels HVG-matched to the real sample.
set -euo pipefail
source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/_env.sh
P=/dcs04/hicks/data/Jan/sim_project/albis_paper
CODE="$P/code/comparison_methods/statistics"
IN="$P/data/figure_5/figure_5B/inputs"
REAL="$P/data/real_data_qc/brain_samples_xenium"
OUT="$P/data/figure_5/figure_5B/real_brain_xenium_hvg_matched"
export MPLCONFIGDIR=/tmp/real_brain_mpl_${SLURM_JOB_ID:-local} OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2

# short labels: they appear inside the JSD box next to the legend
declare -A LABELS=([brain_healthy]="Xenium healthy" [brain_glioblastoma]="Xenium GBM" [brain_alzheimers]="Xenium AD")
for sample in brain_healthy brain_glioblastoma brain_alzheimers; do
  label="${LABELS[$sample]}"
  "$PYTHON_BIN" "$CODE/plot_real_vs_three_methods.py" --real "$REAL/$sample/${sample}_qc.h5ad" --real-label "$label" \
    --albis "$IN/albis_cell_slice4_qc.h5ad" --spider "$IN/spider_cell_slice4_qc.h5ad" \
    --sccube "$IN/sccube_cell_slice4_qc.h5ad" --output-dir "$OUT/$sample"
done
echo "[done] $OUT"
