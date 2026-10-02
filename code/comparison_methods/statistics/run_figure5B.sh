#!/usr/bin/env bash
#SBATCH --job-name=figure5B
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/statistics/logs/figure5B_%j.out
# Figure 5B with Figure 2's code, on Figure 5A data, stored slice 4 (5th from the bottom):
#   1. subset each method/modality to the slice      (subset_slice.py)
#   2. QC with Figure 2's step00_qc_filter.py              (total > 0, >= 3 genes)
#   3. plot_three_methods.py: ALBIS, SPIDER and scCube on the same panels with Figure 2's
#      statistics (count_distribution.py functions), full_panel + qc_filtered. scCube is
#      log-normalized, so it appears on genes detected, sparsity and (cell) log1p only.
# All methods have 556 genes, so Figure 2's panel-matching modes are not needed.
set -euo pipefail
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
P=/dcs04/hicks/data/Jan/sim_project/sim_paper
CODE="$P/code/comparison_methods/statistics"
SRC="$P/data/figure_5/figure_5A_600k"
OUT="$P/data/figure_5/figure_5B"
IN="$OUT/inputs"
SLICE=4
export MPLCONFIGDIR=/tmp/figure5B_mpl_${SLURM_JOB_ID:-local} OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
mkdir -p "$IN"

for modality in cell bin spot; do
  for method in albis spider sccube; do
    raw="$IN/${method}_${modality}_slice${SLICE}.h5ad"
    "$PYTHON_BIN" "$CODE/subset_slice.py" --input "$SRC/$method/$modality.h5ad" --output "$raw" --slice-id $SLICE
    "$PYTHON_BIN" "$P/code/clustering/step00_qc_filter.py" --modality "$modality" --input "$raw" --output "${raw%.h5ad}_qc.h5ad"
  done
  for mode in full_panel qc_filtered; do
    suffix=""; [[ $mode == qc_filtered ]] && suffix="_qc"
    albis="$IN/albis_${modality}_slice${SLICE}${suffix}.h5ad"
    "$PYTHON_BIN" "$CODE/plot_three_methods.py" --modality "$modality" --albis "$albis" \
      --spider "$IN/spider_${modality}_slice${SLICE}${suffix}.h5ad" \
      --sccube "$IN/sccube_${modality}_slice${SLICE}${suffix}.h5ad" --output-dir "$OUT/$modality/$mode"
  done
done
echo "[done] $OUT"
