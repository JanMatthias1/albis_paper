# Supplementary Figure 6: scCube spot aggregation

Code: `sim_paper/code/comparison_methods/overview/figure_supp_6/plot_figure_supp_6.py`.
Outputs: `sim_paper/data/comparison_methods/figure_supp_6/`.

The selected existing figure is `slice5_spot_aggregation.png` (also PDF and SVG).
Its plots, per-spot counts, memberships, aggregation summary and original
provenance were copied unchanged from the archived Figure 5A supplement.
`relocation.json` records that archive location and original file checksums.

The script now writes to the Figure supp 6 output directory. Its default input
is the archived `overview_figure2_600k_35881874` tissue recorded in the existing
figure's provenance, not the newly regenerated Figure 5A tissue. This relocation
does not change the aggregation or display implementation or regenerate data.

Historical reproduction command, from the project root:

```bash
MPLCONFIGDIR=/tmp/figure_supp_6 comparison_methods/env/analysis/bin/python sim_paper/code/comparison_methods/overview/figure_supp_6/plot_figure_supp_6.py
```

The inherited geometry-only analysis displays targets 3, 100, 1000, 5000 and
10000 cells/spot for slice 5 (stored ID 4). Its vectorized grid calculation and
dominant-cell-type summaries are retained as originally implemented. These are
historical derived occupancy summaries; the command does not train a model or
regenerate expression. The saved figure contains 59,699 input cells.
