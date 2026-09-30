# Supplementary Figure 6: scCube spot aggregation

How scCube's spots change with the requested cells per spot, on Figure 5A's
scCube tissue, slice 5 (stored slice_id 4, 59,699 cells).

Code: `plot_figure_supp_6.py`. Outputs: `sim_paper/data/figure_5/figure_supp_6/`.

```bash
sbatch --mem=64G --time=12:00:00 --wrap \
  "comparison_methods/env/sccube/bin/python -u sim_paper/code/comparison_methods/overview/figure_supp_6/plot_figure_supp_6.py"
# redraw from saved tables only: add --redraw-only
```

Native scCube only: for each target (3, 100, 1,000, 5,000, 10,000 cells per
spot) it calls `generate_spot_data_random(platform="Visium", n_cell=target)`
on the Figure 5A scCube cells (scCube's own grid coordinates, as in
`../generate.py`) and `calculate_spot_prop` for per-spot proportions. Figure 5A
uses 3 for bins (ST layout) and 10 for spots (Visium). As a check, target 10
must reproduce the saved Figure 5A slice-5 spots exactly.

Our conventions: spot colour = most abundant type (ties to the lowest type
number, as in Figure 5A); glyph size scales with grid spacing and is not a
physical capture footprint.

Archive: the previous version (a vectorised re-implementation of scCube's grid
on the archived 2,000-gene tissue) is in `../../_archive/figure_supp_6_20260929/`;
its outputs are in `data/figure_5/archive/_archive/figure_supp_6_reimplemented_old_tissue_20260929/`.
