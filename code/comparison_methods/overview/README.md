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

## Native 3D generation and sectioning

Each method builds its own 3D tissue and cuts it into sections with its own
functions; `generate.py` adds no geometry or slicing of its own:

- ALBIS: `simulate_3d_molecule_sphere_multires` (sphere, native sectioning).
- scCube: `generate_pattern_random(spatial_dim=3, is_split=True,
  split_coord="point_z", slice_num=n_slices)` (cube; labels 1-based, stored 0-based).
- SPIDER: `simulate_10X_3d` (cube), then `spider.slice_anndata_by_z` with
  explicit edges `linspace(0, extent_um, n_slices + 1)`; `z_bins=<int>` would
  drop the max-z cell.

Per-slice bins/spots use each method's own 2D capture routine. Bin/spot colour
is the most common type (ties go to the lowest type number, as in ALBIS) of the
method's own composition: scCube `calculate_spot_prop`; SPIDER's returned `W`
for square bins. For SPIDER circular spots `W` is computed but not returned, so
we repeat SPIDER's own one-line product (the single user-approved exception in
`../AGENTS.md`). SPIDER's transition target comes from its
`make_transition_matrix("attractive", 8, strength=0.7)`.

SPIDER is not seedable through `simulate_10X_3d`. Its annealer builds
`AnnealingConfig()` without `random_state`, so cell-type layouts differ between
runs even with the same seed. Coordinates are reproducible, and the neighbour
statistics are stable (8-NN same-type fraction about 0.69 across repeated runs).
The seedable `spider.simulate_cells` was tested and rejected: with default
settings it reaches only 0.19–0.33 against the 0.70 target. Figure 5A
reproducibility therefore rests on the saved SPIDER outputs. The current
600k scCube/SPIDER data predate this change (sections came from
`floor(z / extent * n_slices)`), but both native slicers reproduce those
assignments for all 600,000 cells (0 mismatches, checked 2026-09-28), so they
were not regenerated.

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

## Spot aggregation supplement

[spot_aggregation/](spot_aggregation/) contains the slice5 scCube occupancy-sensitivity figure code and reproduction instructions.

## Intact tissue column

The combined overview now shows Tissue Simulation → Stacked Tissue Slices → Slice5. The first column renders all native continuous3D cell positions; it does not reconstruct tissue from section planes. scCube retains its native-coordinate record and is validated against the uniformly scaled display coordinates. All methods must retain more than10 distinct Z coordinates. Cell types use the same palette in all panels. Previous overview exports are preserved under figure_5A_600k/archive_before_intact_tissue_20260928/.
