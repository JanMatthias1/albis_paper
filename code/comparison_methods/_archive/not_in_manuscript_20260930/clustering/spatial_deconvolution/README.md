# 600k SPIDER deconvolution and three-seed comparisons

Active SPIDER runs use 600,000 cells, ten slices, eight cell types, gene panels
556 and 2,000, and simulation seeds 2025, 101, and 202. Previous pilots and the
superseded three-panel comparison are preserved under `archive/` in the code
and data directories. scCube deconvolution remains archived.

## Selected figure (2026-09-28)

The Figure 5 panel is now `albis_16um_spot_parameters/plot_deconvolution.py`
(ALBIS on the Figure 5A shared tissue, 3 simulation seeds; see that folder's README).
The section below describes the superseded Figure 3-based comparison.

## Recreate the superseded single-axis figure (Figure 3 ALBIS data)

Run from `/dcs04/hicks/data/Jan/sim_project`:

```bash
comparison_methods/env/analysis/bin/python sim_paper/code/comparison_methods/clustering/spatial_deconvolution/plot_seed_mean_comparison.py
```

Reads saved RCTD predictions and truth; no refitting. Writes PNG, PDF, SVG,
per-seed statistics, mean/SD statistics, and provenance into
`data/figure_5/spatial_deconvolution/albis_spider_error_comparison/`.
For each cell type, calculate mean absolute fraction error across scored spots
within each run, then the equally weighted mean and sample SD across three runs.
Faint points show individual runs; all ten slices contribute.

ALBIS uses existing 556-gene results from
`data/figure_3/spatial_deconvolution/RCTD/strong_mix{,_seed101,_seed202}`:
reference seeds vary but the query tissue is fixed. SPIDER uses simulation seeds.
ALBIS truth is molecule fractions; SPIDER truth is contributing-cell fractions.
Type indices are not biologically matched. The comparison is descriptive.

## Active scripts

- `plot_seed_mean_comparison.py`: selected single-axis ALBIS/SPIDER figure.
- `final_spider_figures.py`: SPIDER overview, cell-type errors, and circular-marker
  Slice 5 truth/estimate panels for each seed; outputs `final_spider600k/`.
- `summarize_spider600k.py`: six-run metrics and seed overview; outputs `spider600k_summary/`.
- `prepare_spider600k.py`: native circular SPIDER aggregation of saved 600k cells
  and construction of RCTD input files.
- `run_spider600k.sh`: six-task RCTD job array. Completed job: 35989955.

The other plotting scripts use the same analysis Python environment and accept
no arguments. Run them using the same command prefix as above.

## Inputs and fitting

Sources: `data/figure_5/clustering/cell_type/experiment_600k/runs/`
`spider_genes{556,2000}_seed{2025,101,202}/`.
Outputs: `data/figure_5/spatial_deconvolution/`
`spider600k_genes{556,2000}_seed{2025,101,202}/`, including inputs, logs,
input audits, fitted RCTD objects, estimates, and metrics.

Native circular spots use diameter 55 µm and spacing 100 µm across all ten slices.
The original 10k-cell Splatter pool supplies the reference; it contributed to
query generation and is not held out. No rounding, normalization, or expression
rescaling is added. Preparation imports metadata helpers from
`code/comparison_methods/overview/generate.py` and calls native SPIDER aggregation.
Fitting reuses `code/clustering/spatial_deconvolution/weak_domain/run_rctd_spot.R`
with full mode, default QC, and fixed R seed 0. Environments: project
`comparison_methods/env/spider` for preparation and `sim_paper/env/rctd` for fitting.

For a new fit, submit `sbatch --array=0-5%2` followed by the full path to
`run_spider600k.sh`, after source generation completes. Existing input directories
are deliberately protected against overwrite; current completed runs need no rerun.
