# Figure 5A: 600,000-cell 3D overview

Output: `sim_paper/data/figure_5/figure_5A_600k/`. Everything in that
folder is regenerated there by one command; nothing is linked from older runs.

## Reproduce

From anywhere (submits five SLURM jobs; the output folder must not exist):

```bash
bash sim_paper/code/comparison_methods/overview/submit_figure5A.sh [OUT]
```

| Stage | What | Environment |
|---|---|---|
| splatter | 556-gene Splatter expression pool, seed 20260921, 10,000 cells | `comparison_methods/env/splatter` |
| albis | ALBIS shared tissue: `applications_albis/cross_modality_alignment/strong_domain_mix/generate_strongmix_offsets.py` with its defaults (Figure 4C strong-mix 16 µm parameters, seed 2025, 600k cells, 556 genes) | `sim_paper/env/albis-tutorial` |
| spider | `generate.py --method spider` | `comparison_methods/env/spider` |
| sccube | `generate.py --method sccube` | `comparison_methods/env/sccube` |
| finish | `verify_reproducibility.py`, then `render.sh` | `comparison_methods/env/analysis` |

Redraw only (no simulation): `bash sim_paper/code/comparison_methods/overview/render.sh [OUT]`.

## Files

- `submit_figure5A.sh`: the full pipeline above.
- `settings_figure5A.json`: all Figure 5A settings (copied into OUT as `settings.json`).
- `generate.py`: scCube or SPIDER, one native 3D tissue then native slices, bins
  and spots. Also used by `../native556/` and `../clustering/cell_type/run_task.py`.
- `verify_reproducibility.py`: compares the new Splatter pool with
  `archive/overview_native556_20260926/splatter` and the new ALBIS outputs with Figure 4
  `strong_domain_mix_shift3x/data`; writes `OUT/reproducibility_check.json`.
- `plot_combined.py`, `render.sh`: the figure.
- `figure_supp_6/`: supplementary Figure 6, scCube spot aggregation vs cells per
  spot on Figure 5A slice 5, using scCube's own `generate_spot_data_random` and
  `calculate_spot_prop`. Outputs: `sim_paper/data/figure_5/figure_supp_6/`.

## What is native

All three methods build the whole 3D tissue before slicing (ALBIS a sphere,
scCube and SPIDER a cube). Slicing only labels cells; bins and spots are then
made per slice. ALBIS slices molecules; scCube and SPIDER slice whole cells.

| Step | scCube | SPIDER |
|---|---|---|
| Expression | `train_vae_and_generate_cell` (from the Splatter pool) | `get_sim_cell_level_expr` (from the Splatter pool) |
| 3D layout | `generate_pattern_random(spatial_dim=3)` | `simulate_10X_3d` |
| Neighbour target | – | `make_transition_matrix("attractive", 8, strength=0.7)` |
| Slicing | `is_split=True, split_coord="point_z"` (1-based, stored 0-based) | `slice_anndata_by_z` with edges `linspace(0, extent, n+1)` (`z_bins=<int>` drops the max-z cell) |
| Bins / spots | `generate_spot_data_random` (ST / Visium) | `sim_expr.get_sim_spot_level_expr` (square / circle) |
| Per-capture composition | `calculate_spot_prop` → `obsm["sccube_spot_prop"]` | returned `W` → `obsm["spider_W"]`; circle spots: see exception |

Our conventions (not simulation): bin/spot colour is the most common type of the
method's own composition (ties go to the lowest type number, as in ALBIS);
scCube grid units are scaled by `extent_um / sccube_grid_size` (512.5 µm);
bin/spot z is the slice midpoint (unused by any calculation or plot); the
Splatter pool is the expression input for scCube and SPIDER.

Exception (user-approved, `../AGENTS.md`): SPIDER's circle branch computes `W`
but returns only three values, so we repeat its own line
`spot_cell_idx_matrix * onehot_ct` with SPIDER's membership matrix and
`spider.core.get_onehot_ct`.

## Reproducibility

- Splatter and ALBIS are deterministic; `reproducibility_check.json` records
  byte/array identity with the earlier runs.
- scCube is seeded (`set_seed`, torch seed).
- SPIDER cannot be seeded through `simulate_10X_3d`: its annealer builds
  `AnnealingConfig()` without `random_state`. Coordinates reproduce; cell-type
  layouts differ between runs, with stable statistics (8-NN same-type ≈ 0.69).
  The seedable `spider.simulate_cells` reaches only 0.19–0.33 of the 0.70 target
  with default settings and was rejected. The saved outputs are the reference.

## Genes

All methods have 556 genes. ALBIS's are its own marker/shared/noise genes;
scCube/SPIDER reproduce the 556 Splatter genes. Counts match, identities do not.

## Archive

- Code: `../_archive/overview_20260929/` (superseded submit scripts, old logs,
  `generate.py` with the removed ALBIS/`--method all` path; `moved_files.json`).
  Older: `../_archive/overview_20260927/`.
- Data: `data/figure_5/archive/_archive/figure_5A_superseded_20260929/`
  (2,000-gene figure, 2026-09-28 native 2,000-gene run, old pools and trial runs).
