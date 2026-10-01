# Comparison with other spatial simulators (Figure 5)

Compares ALBIS with two 3D spatial simulators, [scCube](https://github.com/ZJUFanLab/scCube)
and [SPIDER](https://github.com/YANG-ERA/Spider), on the same 600,000-cell tissue
size and 556 genes. Every method runs through its own public API; outputs are
saved as each method produces them. Outputs go to `data/figure_5/`.

## 5A: 3D overview (`overview/`)

Each method builds one 3D tissue (ALBIS a sphere, scCube and SPIDER a cube),
which is cut into 10 slices and captured as cells, bins and spots. scCube and
SPIDER take their expression from a shared Splatter reference; ALBIS uses its
own genes.

- `submit_figure5A.sh`: full pipeline (Splatter → ALBIS → SPIDER → scCube → checks → figure)
- `generate.py`: native scCube / SPIDER generation, slicing, bins and spots
- `verify_reproducibility.py`: checks that Splatter and ALBIS outputs match the earlier runs exactly
- `plot_combined.py`, `render.sh`: the figure (`render.sh` redraws without simulating)
- `settings_figure5A.json`: all settings

## 5B: expression statistics (`statistics/`)

The Figure 2 count-distribution comparison (`code/count_distribution/`),
applied to slice 5 of the 5A data: ALBIS vs SPIDER and ALBIS vs scCube. scCube
outputs log-normalized values rather than counts, so `sccube_compare.py`
compares only what applies (genes detected, sparsity, normalized values).

- `run_figure5B.sh`: subset to one slice → QC → comparisons

## 5D: capture close-up (`figure_5D/`)

The central 290 × 290 µm of slice 5 for each method, with a common hexagonal
spot overlay: ALBIS molecules, scCube and SPIDER cells. SPIDER uses a separate
seed here, because with the 5A seed its cell positions coincide with scCube's.

- `submit_distinct_seeds.sh`: Splatter check → SPIDER run → `actual_capture.py`
- `native_spots.py`: each method's own spots; `plot_albis_3d.py`: 3D view of the ALBIS molecules

## Supplementary Figure 6: scCube spot aggregation (`overview/figure_supp_6/`)

How scCube's spots change with the requested cells per spot (3 to 10,000) on
the 5A scCube tissue, using scCube's own `generate_spot_data_random`.

## Compute: runtime and memory (`compute/`)

Runtime and memory of ALBIS, scCube and SPIDER from 10k to 5M cells, on one
CPU thread, for three seeds. Setup (Splatter reference, scCube VAE training)
is timed separately from generation.

- `submit_seeds.sh` → `prepare_seeds.py`, then `native_extension.sbatch` (ALBIS, SPIDER) and `sccube_seed_batch.sbatch` (scCube)
- `albis_native.py`, `spider_native.py`, `sccube_split.py`: timed native calls
- `report_large.py`, `plot_compute_four_panel_preliminary.py`: tables and figures
- `METHODS.md`: methods text

## Conventions

These choices are ours, not outputs of the methods:

- Bins and spots are coloured by the most common cell type in each method's own
  composition (ties go to the lowest type number).
- scCube grid units are converted to µm (512.5 µm per unit); native values are kept.
- SPIDER's circular-spot branch computes the composition matrix `W` but does not
  return it, so `generate.py` repeats SPIDER's own line to recover it (used only
  for 5A colours).
- SPIDER's cell-type layout cannot be seeded through `simulate_10X_3d`, so its
  layouts differ between runs; positions and summary statistics reproduce.
- Single-slice panels use stored `slice_id` 4, shown as "Slice 5".

## Environments

Separate conda environments for Splatter (R), scCube, SPIDER (`st-spider`)
and the comparison plots; `env/albis-tutorial` for ALBIS.
