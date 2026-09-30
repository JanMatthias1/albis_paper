# Figure 5B: expression-distribution comparison

Run `sbatch run_figure5B.sh` to regenerate all three modalities and both gene-panel variants from the current `data/comparison_methods/figure_5A_600k/{albis,spider,sccube}/` files. No simulation is regenerated or modified. Uses displayed Slice 5 (stored slice_id=4), matching Figure 5A.

Outputs: `data/comparison_methods/figure_5B/{cell,bin,spot}/distribution_comparison{,_qc_hvg}/`. Each contains mean–variance, mean–dropout, sparsity, total expression, detected genes, raw/normalized/log-expression plots, a legend, and comparison_summary.json with input paths and numerical summaries.

Both variants exclude captures explicitly marked empty or unassigned, preserving the original script's filtering behavior. Native panels have 556 ALBIS genes versus 2,000 competitor genes. The `_qc_hvg` variant retains 556 genes per method, selecting competitors' genes by variance after normalize_total/log1p while computing distribution statistics from the original expression values. This matches panel size, not gene identity. scCube continuous output is not rounded. Differences in expression/noise models and capture counts remain.

The former Figure 4A outputs are preserved under `data/comparison_methods/archive/figure_4A_pre_figure5B_20260927/`; the old script is under `archive/pre_figure5B_20260927/`. The new plot code still imports the shared `code/count_distribution/count_distribution.py` statistics helpers.

Initial regeneration job: 35986757.

## Log-normalized cell mean–variance plot

`plot_log_normalized.py` adds `mean_variance_log_normalized.{png,pdf,svg}`, a per-gene CSV and transformation/selection metadata JSON to both cell output directories. Run independently with `sbatch run_log_normalized.sh`; the full run_figure5B.sh also includes it.

ALBIS/SPIDER are normalized to 10,000 using their full native panels and log1p-transformed. scCube's native reconstructed log-normalized values are used unchanged. Linear plot axes show means/variances of these already transformed values, without Poisson or NB curves. The matched version selects the 556 highest-variance genes on this comparison scale, without renormalizing after selection. Its scCube selection therefore differs from the historical raw-panel selection, which transformed scCube a second time. Native panel sizes and gene identities remain different. Bin/spot panels are not transformed by this script.

Cell-panel correction: `raw_norm_log_compare.png` now shows ALBIS/SPIDER in raw and normalized panels, and all three methods in log1p(normalized). scCube values are unchanged and are never normalized/logged again. Cell-level matched-panel selection now also uses scCube's native variance. `rerun_cell.sh` regenerates the cell plots and summaries affected by this corrected selection. This change is cell-only; aggregate panels remain separate legacy diagnostics.

## Corrected aggregate scale handling

Across all modalities, scCube is never normalized/logged again, including during variance-based panel selection. Count-based mean–variance/NB, mean–dropout and total-count panels contain ALBIS/SPIDER only. `sccube_native_expression.png` describes native scCube mean–variance and totals separately. Sparsity and nonzero-gene panels remain descriptive three-method comparisons; exact zeros depend on the output representation.

For cells, scCube appears unchanged in the log1p(normalized) histogram. For bins/spots it appears unchanged in a fourth panel labelled summed cell log-expression; it is not compared by JSD against count-derived transformations. All pairwise annotations use `ALBIS vs METHOD`. NB dispersion for scCube is null in summaries. Historical field names `genes_per_cell_*` refer to observations at the selected modality; `native_total_expression_median` records all methods, while `total_counts_*` now excludes scCube.

## JSD overview across ten slices

`jsd_overview.py` computes one QC/HVG log-expression JSD per slice for ALBIS vs SPIDER (cell/bin/spot), and ALBIS vs scCube (cell only). `run_jsd_overview.sh` runs a slice array, followed by the same wrapper with `aggregate`. Initial jobs: 35987420 (array), 35987483 (dependent aggregation).

Outputs: `data/comparison_methods/figure_5B/jsd_overview_qc_hvg/jsd_boxplot.{png,pdf,svg}`, `jsd_by_slice.csv`, `protocol.json`, and per-slice samples/metadata. Boxes summarize ten slices from one tissue realization, not independent simulation seeds. Shared 60-bin histograms and reproducible positive-entry sampling are used for all comparisons. Current QC/HVG composite selection/transformation order is retained. The plot is an additional overview; existing panels are retained.

## Current ALBIS source

Figure 5B uses exactly the Figure 5A source links: Figure 4
`cross_modality_alignment/strong_domain_mix_shift3x/data`, with
`cell/simulation_cell_z.h5ad`, `bin16um/simulation_bin_z.h5ad`, and
`spot/simulation_spot_z.h5ad`. All three arise from the shared 600k-cell tissue.
The previous Figure 2-calibrated plots are preserved in the data archive;
`figure_5B/source_change.json` records their exact location. This source update
retains all existing filtering, gene selection and scale-handling rules.
