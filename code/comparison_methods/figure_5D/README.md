# Figure 5D: actual simulator coordinates under a hexagonal overlay

Run (three SLURM jobs; the output folder must not exist):

```bash
bash sim_paper/code/comparison_methods/figure_5D/submit_distinct_seeds.sh [OUT]
```

Default output: `sim_paper/data/figure_5/figure_5D_actual_coordinates/distinct_seeds/`.

| Stage | What |
|---|---|
| splatter | 556-gene pool, seed 20260921; the job fails unless it is byte-identical to Figure 5A's |
| spider | `overview/generate.py --method spider`, Figure 5A settings with seed 20260922, cells only |
| plot | `actual_capture.py` (defaults = this output folder) |

Redraw only: `sim_paper/env/albis-tutorial/bin/python actual_capture.py`.

## What is shown

The central 290 × 290 µm field of slice 5 (stored `slice_id` 4) with a common
seven-spot hexagonal overlay (55 µm diameter, 100 µm spacing). The overlay is
illustrative; geometric inclusion is not any method's native spot membership.
Counts are molecules for ALBIS and cells for scCube/SPIDER, not comparable units.

- ALBIS: native molecules regenerated from the exact Figure 5A config and seed
  (source hash and all 600k cell coordinates/labels checked against the saved
  data); cached in `figure_5D_actual_coordinates/albis_molecules_roi.npz`.
- scCube: saved Figure 5A cell positions.
- SPIDER: cell positions from the separate seed-20260922 run. With Figure 5A's
  shared seed, SPIDER's positions equal scCube's (see `../README.md`).

`plot_albis_3d.py`: 3D view of the cached ALBIS molecules and cells
(`figure_5D_actual_coordinates/figure5d_albis_3d.*`); no regeneration.

## Archive

`../_archive/not_in_manuscript_20260930/figure_5D/`: the synthetic schematic
(`plot_hexagonal_capture.py`) and the shared-position runner
(`run_actual_capture.sh`). Their outputs are in
`data/figure_5/archive/not_in_manuscript_20260930/figure_5D_hexagonal_capture/`
and `figure_5D_actual_coordinates/archive_shared_positions_20260930/`.
