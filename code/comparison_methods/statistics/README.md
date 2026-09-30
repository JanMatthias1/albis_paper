# Figure 5B: expression-distribution comparison

Figure 2's method, applied to the Figure 5A data: compare distributions and
compute JSD. Nothing is simulated or changed.

```bash
sbatch sim_paper/code/comparison_methods/statistics/run_figure5B.sh
```

Inputs: `data/figure_5/figure_5A_600k/{albis,spider,sccube}/{cell,bin,spot}.h5ad`,
stored slice 4 (the 5th of 10 slices from the bottom, shown as "Slice 5" in
Figure 5A; Figure 2 uses the same slice since 2026-09-29). Environment: Figure 2's
(`code/count_distribution/_env.sh`).

Steps (`run_figure5B.sh`):

1. `subset_slice.py`: keep slice 4 of each file (Figure 2's script restricts only
   its first dataset to a slice, so both sides are pre-subset).
2. QC with Figure 2's `code/clustering/step00_qc_filter.py` (total > 0, ≥ 3 genes).
3. **ALBIS vs SPIDER**: Figure 2's `count_distribution/count_distribution.py`,
   unchanged, in its `full_panel` (no QC) and `qc_filtered` modes. All methods
   have 556 genes, so the panel-matching modes are not needed. Only
   `--title-context "ALBIS vs SPIDER"` changes the plot titles.
4. **ALBIS vs scCube**: `sccube_compare.py`, built from Figure 2's own functions
   (`compute_jsd`, bins, sampling). scCube outputs log-normalized values, not
   counts, so count panels do not apply: genes detected per observation and
   sparsity for all modalities; the log1p(normalized) histogram with JSD for
   cells only (ALBIS normalized as in Figure 2; scCube unchanged); for bins and
   spots (sums of cell log-expression) a descriptive native panel, no JSD.

Outputs: `data/figure_5/figure_5B/<modality>/albis_vs_{spider,sccube}/{full_panel,qc_filtered}/`
plus `figure_5B/inputs/` (slice-4 and QC'd inputs). As in Figure 2, read
mean-variance/dropout from `full_panel`; for bins and spots use `qc_filtered`.

Archive: previous 5B code (own plotting copies, custom variance matching, JSD
across 10 slices, log-normalized mean-variance) in
`archive/pre_figure2_code_20260929/`; its outputs in
`data/figure_5/archive/figure_5B_before_figure2_code_20260929/`.
