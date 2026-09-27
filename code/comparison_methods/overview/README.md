# Figure 5A — 600,000-cell 3D overview

Current outputs and configuration:
`sim_paper/data/comparison_methods/figure_5A_600k/`.

## Redraw the current figure

From the project root:

```bash
bash sim_paper/code/comparison_methods/overview/render.sh
```

This reads the saved inputs and settings and rewrites
`overview_combined.{png,pdf,svg}` plus the plotting summary. It does not run
simulations. Titles are “3D Tissue Simulation” and “Slice 5” (stored ID 4),
with a cell-type legend and no footer.

## Active files

- `plot_combined.py`: the current figure layout, colors and rendering.
- `render.sh`: redraws the canonical figure; optional first argument selects
  another prepared input directory.
- `generate.py`: native simulation worker, including st-spider 1.2.0 and
  aggregation validation. Also used by the separate native556 workflows.
- `submit_spider600k.sh`: cluster generation and rendering for a **new,
  prepared** directory containing settings, ALBIS/scCube links and a Splatter
  reference. It refuses to overwrite an existing SPIDER dataset.
- `logs/`: the completed 600k SPIDER job log.
- `submit_sccube600k.sh`: prepares a new directory, trains scCube and generates
  600,000 cells plus bins/spots, then renders the figure. Reuses the current
  synthetic Splatter pool and links the current ALBIS/SPIDER datasets.

To regenerate scCube without overwriting the current figure:

```bash
sbatch sim_paper/code/comparison_methods/overview/submit_sccube600k.sh \
  /dcs04/hicks/data/Jan/sim_project/sim_paper/data/comparison_methods/figure_5A_sccube_rerun
```

The output directory must not exist. This trains the VAE afresh using the
saved 600k settings; it does not regenerate the reference pool or other methods.

## Inputs and dependencies

All methods have 600,000 cells and ten sections. ALBIS uses the uncropped
`figure_4/cross_modality_alignment/strong_domain_mix_shift3x/data` files,
with aligned coordinates and 16 µm bins. scCube reuses the completed 600k
realization. SPIDER was generated with installed `st-spider==1.2.0`.
Exact paths and parameters are in the output directory's `settings.json`,
`input_provenance.json` and method manifests. Do not move linked source data.

Plotting imports `code/manuscript_style.py` for default colors; the current
stronger palette, opacity and spot sizes are specified in `settings.json`.
Python environments live in `comparison_methods/env/`. Full generation also
uses shared helpers and Splatter code under `comparison_methods/code/` and
the local ALBIS package; this folder is not a standalone software bundle.

## Cleanup

Superseded preparation scripts, renderers, settings, launchers and the old
README are archived in `../_archive/overview_20260927/`, with a move manifest.
Archived scripts preserve historical code and may reference their original
locations; they are not current entry points.
The 556-gene experiment launchers now live in `../native556/` and statistical
distribution plotting in `../statistics/`. These serve other comparison panels.
