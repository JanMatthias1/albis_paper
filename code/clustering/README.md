# Clustering workflow (Figure 3)

Two conda envs are involved:

```bash
# Everything except the BANKSY build step
conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial
python -m pip install -r sim_paper/code/clustering/requirements.txt

# BANKSY build step only (01_build_banksy_matrix.py) -- banksy_py pins an
# older scanpy/numpy/anndata/pandas/scikit-learn/scipy stack, incompatible
# with albis-tutorial. See sim_paper/env/create_banksy_env.sh.
conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/sim-app-banksy
```

All scripts here source `_env.sh` (albis-tutorial) or `_env_banksy.sh`
(sim-app-banksy) for the conda-activation boilerplate.

## How to (re)build each `data/figure_3/` output

Quick reference: which script(s) to run, in order, for each data product.

**`pca_harmony_single_cell/`** -- one self-contained script per modality,
covered in the next section:
```bash
sbatch sim_paper/code/clustering/pca_harmony_single_cell/cell_celltype_panel.sh
sbatch sim_paper/code/clustering/pca_harmony_single_cell/bin_celltype_panel.sh
sbatch sim_paper/code/clustering/pca_harmony_single_cell/bin16um_celltype_panel.sh
sbatch sim_paper/code/clustering/pca_harmony_single_cell/spot_celltype_panel.sh
```

**`banksy_batch_compare/`** (data dir; the code lives at `misc/banksy_batch_compare/`,
not its own top-level folder) -- build, then patch spot with corrected
dispersion (depends on `cellbin_batch_sigma_slide/` below having already
built spot's corrected-dispersion data), then add the qualitative plots:
```bash
sbatch sim_paper/code/clustering/misc/banksy_batch_compare/banksy_batch_compare.sh
# then, once cellbin_batch_sigma_slide/spot_batch_sigma_slide_corrected.sh (below) has landed --
# overwrites spot/{prebatch,tuned} only, current canonical spot data:
sbatch sim_paper/code/clustering/misc/banksy_batch_compare/rebuild_spot_batch_compare_corrected_params.sh
sbatch sim_paper/code/clustering/misc/banksy_batch_compare/run_banksy_batch_compare_umap.sh
sbatch sim_paper/code/clustering/misc/banksy_batch_compare/run_banksy_batch_compare_true_vs_pred.sh
```
Note bin's official prebatch/tuned pair does NOT live here (this runs bin at
the stale k_geom=200) -- it's built by `cellbin_batch_sigma_slide/` instead,
see below and `plot_banksy_batch_compare_3mod.py`'s header.

**`cellbin_batch_sigma_slide/`** -- 2026-09-16: pared down to just the
seed-replication + plotting + UMAP/Harmony-diagnostic code (by decision --
the original lambda/k_geom-tuning and base-grid generation scripts that
first built this directory moved to `misc/`, see below; the data they built
is untouched on disk, this folder just no longer carries the from-scratch
recipe for it). What's here now:
```bash
# seed-replication (2 extra seeds/point across all 3 modalities + 2 new spot points beyond canonical):
sbatch sim_paper/code/clustering/cellbin_batch_sigma_slide/batch_sigma_slide_seeds.sh
# UMAP / before-after-Harmony diagnostics:
sbatch sim_paper/code/clustering/cellbin_batch_sigma_slide/run_batch_sigma_slide_umap.sh
sbatch sim_paper/code/clustering/cellbin_batch_sigma_slide/cell_umap_gap.sh
```
`gen_batch_sigma_slide_seed_tasks.py` (deterministic task-list generator) and
`batch_sigma_slide_seed_tasks.tsv` (its output) sit alongside
`batch_sigma_slide_seeds.sh` and feed it directly.

If you need to rebuild the base grid (not the seed replicates) from scratch --
e.g. the data was deleted -- the original generation chain is now in `misc/`:
`cellbin_batch_sigma_slide.sh` -> `bin16um_kgeom100_endpoints.sh` (depends on
`banksy_batch_compare.sh`'s QC'd bin h5ad above) -> `batch_sigma_slide_fine.sh`
-> `spot_batch_sigma_slide_corrected.sh` (current canonical spot/, run last --
it rebuilds all 16 spot points and supersedes the spot ones the earlier steps
built). `misc/spot_batch_sigma_slide_weakmix.sh` (-> `spot_weakmix/`) is a
separate validation sibling, not required for the main curve either way.
Gotcha either way: cell's `bs0`/`bs1.5` and spot's `bs0`/`bs0.3` canonical
points are NOT built by any of this -- `plot_batch_sigma_slide_final.py`
pulls those from `banksy_batch_compare/` instead (build that dir first).

**These two directories cross-depend on each other** (see the dependency
notes above) -- neither is fully reproducible standalone.

Then, the manuscript plot itself:
```bash
python sim_paper/code/clustering/cellbin_batch_sigma_slide/plot_batch_sigma_slide_final.py
```

**`ari_recovery_summary/`** -- pure aggregator, no generation step; just
discovers whatever's already on disk under `pca_harmony_single_cell/` (for
`cell_type_true`) and `data/clustering_<tag>/` (for `domain_true` -- not one
of the four dirs above, see "Output layout" below):
```bash
python sim_paper/code/clustering/ari_recovery_summary/plot_ari_recovery.py
```

## Reproduce Figure 3: run one panel script per modality x task

**Cell-type**: the fastest path is one self-contained script per modality,
under `pca_harmony_single_cell/`:

```bash
sbatch sim_paper/code/clustering/pca_harmony_single_cell/spot_celltype_panel.sh
sbatch sim_paper/code/clustering/pca_harmony_single_cell/bin_celltype_panel.sh
sbatch sim_paper/code/clustering/pca_harmony_single_cell/bin16um_celltype_panel.sh
sbatch sim_paper/code/clustering/pca_harmony_single_cell/cell_celltype_panel.sh
```

These read the exact same sim `.h5ad` files as Figure 2 (`data/figure_2/<tag>/`,
tags carrying the per-modality `batch_sigma`: cell 1.5, bin8 0.8, bin16 0.7,
spot 0.3).

**Domain**: NOT one self-contained script per modality in current practice --
see "How to (re)build each `data/figure_3/` output" above for the actual
`banksy_batch_compare/` + `cellbin_batch_sigma_slide/` chain that feeds the
manuscript panel. `legacy/spot_domain_panel.sh`, `legacy/bin_domain_panel.sh`,
`legacy/cell_domain_panel.sh` are an older, self-contained-per-modality
attempt at the same thing (generate its own strong-`domain_type_mix` dataset
-> BANKSY -> ARI -> Leiden, writing to `data/clustering_<tag>/`, not
`data/figure_3/`) -- still runnable, but superseded in practice: the current
manuscript panel (`plot_domain_vs_celltype_recovery.py`) does not read from
their output. Moved to `legacy/` 2026-09-16 to make that clear.

Each is idempotent end-to-end: generates the sim dataset if it isn't already
on disk, QC-filters it, clusters it (plain PCA+Harmony for the cell-type
panels, BANKSY+PCA+Harmony for the domain panels -- see "Why two pipelines"
below), computes resolution-matched ARI against ground truth, and produces
the qualitative ground-truth-vs-predicted UMAP + contingency-heatmap plots.
Rerunning a script after its output already exists skips straight to
whichever step is missing -- but note this checks for the *final* artifact
of each step, not every plot inside it: if you change a step's plotting code
without changing its output data, rerunning the panel script won't
regenerate that step's plots (its h5ad/json output already exists, so the
whole step is skipped). Force a re-plot with the underlying script's
`--plots-only` flag (`pca_harmony.py`, `03_clustering_plots.py`)
against the existing output instead.

Each script's header documents its exact sim config (`batch_sigma`,
weak-vs-strong `domain_type_mix`, packing fraction) and why -- see the script
itself, and figure.md, for the reasoning.

### Why two pipelines

Cell-type recovery and domain recovery are split by design, each on the
dataset built for its own question:

- **Cell-type panels** (`*_celltype_panel.sh`): plain PCA -> Harmony ->
  Leiden, on the weak/manuscript-baseline `domain_type_mix`. This already
  recovers `cell_type_true` well once `batch_sigma` is tuned per modality.
- **Domain panels** (`legacy/*_domain_panel.sh`, superseded -- see above):
  BANKSY -> PCA -> Harmony ->
  Leiden, on the strengthened `domain_type_mix` (`STRONG_DOMAIN_TYPE_MIX`,
  needed because the manuscript-baseline mix makes several domains
  compositionally indistinguishable -- see figure.md, 2026-08-25). BANKSY is
  never used for cell-typing; plain PCA+Harmony is never used for domain
  recovery on the strong-mix dataset (Harmony's slice-batch-correction
  collapses cell-type recovery when combined with the strong mix).

## Steps, if you need to run them individually

1. **Generate** the sim dataset: `sim_paper/code/data/generate_simulation_noisy.py`
   (not in this directory -- see `sim_paper/code/data/`).
2. **`00_qc_filter.py`** -- drops fully-empty and near-empty observations
   before clustering (targets neighbor-graph fragmentation from empty
   bins/spots, not sim-vs-real fairness; see the script's docstring).
   ```bash
   python sim_paper/code/clustering/00_qc_filter.py --modality bin --packing-tag <tag>
   ```
3. **Cluster** -- one of:
   - `pca_harmony.py`: plain PCA -> Harmony (on `obs["slice_id"]`, the batch
     key albis's synthetic batch effects are applied per). Writes PCA and
     before/after-Harmony UMAP diagnostic plots plus PC-pairs plots
     (`pc_pairs.py`, imported, not a standalone CLI).
     ```bash
     python sim_paper/code/clustering/pca_harmony.py --modality bin --input <qc.h5ad> --output <out.h5ad>
     ```
   - `01_build_banksy_matrix.py` (needs the `sim-app-banksy` env): BANKSY ->
     PCA -> Harmony, staggering each slice's coordinates first (all slices
     share one local x/y range, so an unstaggered neighbor graph would treat
     different slices as spatial neighbors). Writes `obsm["X_pca_harmony"]`
     in the same shape `pca_harmony.py` does, so `03_clustering_plots.py`
     needs zero changes to consume either.
     ```bash
     python sim_paper/code/clustering/01_build_banksy_matrix.py --modality cell --packing-tag <tag>
     ```
4. **`03_clustering_plots.py`** -- Leiden or Louvain clustering
   (`--algorithm`, default `leiden`) + UMAP from the pca_harmony pipeline's
   output (`--pipeline` kept as a single-choice flag, `pca_harmony`, since a
   second gene-space Harmony pipeline lived here before being deleted
   2026-08-27; its dead branch was stripped 2026-09-16), plots UMAP
   colored by `domain_true`/`cell_type_true`/`cluster_label`, and (if both a
   `--cluster-key` and a `*_true` ground-truth column are present) a
   ground-truth-vs-predicted UMAP + contingency heatmap.
   ```bash
   python sim_paper/code/clustering/03_clustering_plots.py --modality bin --resolution 0.5
   ```
5. **`02_leiden_resolution_sweep.py`** -- Figure 3B: binary-searches the Leiden
   resolution until the cluster count matches the true category count (so
   ARI isn't confounded by over/under-clustering), separately for
   `domain_true` and `cell_type_true`.
   ```bash
   python sim_paper/code/clustering/02_leiden_resolution_sweep.py --modality bin --packing-tag <tag>
   ```
6. **`plot_ari_recovery.py`** -- grouped bar chart of `cell_type_true` ARI
   across modality, from `pca_harmony_single_cell/`. Trimmed 2026-09-16 to
   drop the `domain_true`/BANKSY half (it only ever read
   `data/clustering_<tag>/<modality>/banksy_ari_recovery/`, populated
   exclusively by `legacy/*_domain_panel.sh`, last built 2026-08-25) -- for
   domain-vs-celltype comparison see `misc/plot_domain_vs_celltype_recovery.py`
   instead, which reads current data.
   ```bash
   python sim_paper/code/clustering/ari_recovery_summary/plot_ari_recovery.py
   ```

Corresponding `run_*.sh` SLURM wrappers for steps 2-6 (`0-2` array over
spot/bin/cell where relevant) live in `misc/` (`run_qc_filter.sh`,
`run_pca_harmony.sh`, `run_clustering_leiden.sh`, `run_ari_vs_ground_truth.sh`)
-- moved from the top level 2026-09-16: nothing in the current pipeline calls
them (every panel/sweep script invokes the `.py` steps directly), and their
only `figure.md` references are old one-off historical reproduce recipes
(Figure 2 write-up, the 2026-09-06 BANKSY tuning sweep), not live
instructions. Still useful if you want to hand-run a custom one-off sweep the
way that work did -- see figure.md's "Reproduce" code blocks for worked
examples.

## Output layout

**Cell-type panels** write under `data/figure_3/pca_harmony_single_cell/<modality>/`:

```
data/figure_3/pca_harmony_single_cell/
  <spot|bin|cell>/
    pca_harmony_qc/
      simulation_<modality>_z_pca_harmony_qc.h5ad
      plots/   # pca_before/after_harmony_by_*.png, umap_pca_harmony_before_after_by_*.png, pc_pairs/
    ari_recovery_qc/
      ari_summary_<modality>.json
    leiden_pca_qc_celltype_matched/
      simulation_<modality>_z_leiden_pca_res<tag>.h5ad
      plots/   # umap_by_*.png, umap_true_vs_predicted_*.png, contingency_*.png
```

**`legacy/*_domain_panel.sh` output is not under `data/figure_3/`** and is
superseded in practice (see "Reproduce Figure 3" above for the actual
manuscript-panel data path, `banksy_batch_compare/` + `cellbin_batch_sigma_slide/`).
Kept for reference: each modality's domain-panel dataset carries its own tag
(differs because `batch_sigma`/`domain_type_mix` differ per modality) and
lives at `data/clustering_<tag>/<modality>/`, mirroring the same three-step
layout (`banksy_pca_harmony_qc/`, `banksy_ari_recovery/`,
`leiden_pca_banksy_domain_matched/`) under a different root:

| modality | domain-panel tag |
|---|---|
| bin | `packing_pf0p04_bsigma05_strongdomainmix` |
| cell | `log_mu_-2.5_theta_0.25_strongdomainmix` |
| spot | `packing_pf0p04_bsigma015_strongdomainmix` |

Check the modality's `legacy/*_domain_panel.sh` header for its current
`SIM_TAG` -- these predate the retune, so treat the table above as a
snapshot, not a guarantee.

`banksy_batch_compare/` and `cellbin_batch_sigma_slide/` -- the data sources
the current manuscript domain panel actually uses -- live correctly under
`data/figure_3/`, per "How to (re)build each `data/figure_3/` output" above.

## `misc/`

One-off exploration/tuning scripts, kept for provenance, not wired into any
canonical panel script. As of 2026-09-16 the scripts that actually produce
`data/figure_3/{pca_harmony_single_cell,ari_recovery_summary,cellbin_batch_sigma_slide}/`
moved out to their own top-level folders (see "How to (re)build each
`data/figure_3/` output" above); `banksy_batch_compare/` code stayed under
`misc/banksy_batch_compare/` (not its own top-level folder -- see that
section above for why and its commands). What's left flat in `misc/`:

- `batch_sigma_sweep.sh`, `batch_sigma_opt_round2.sh`, `tune_banksy_bin_kgeom.sh`,
  `tune_banksy_cell.sh`, `replot_celltype_panels_v2.sh`, `replot_pca_harmony_qc.sh`
  -- 2026-08-27 initial `batch_sigma`/BANKSY tuning, superseded by later work above.
- `banksy_lambda_kgeom_sweep_16um.sh` + `summary_banksy_lambda_kgeom.py` --
  2026-09-08/09 lambda x k_geom re-tune (winner: bin16um k_geom=100, feeds
  `cellbin_batch_sigma_slide/`).
- `cellbin_batch_sigma_slide.sh`, `bin16um_kgeom100_endpoints.sh`,
  `batch_sigma_slide_fine.sh`, `spot_batch_sigma_slide_corrected.sh`,
  `spot_batch_sigma_slide_weakmix.sh` -- moved from `cellbin_batch_sigma_slide/`
  2026-09-16 by decision: the base-grid generation chain that originally
  built the batch_sigma-slide data, kept here for from-scratch reproducibility
  but not needed day-to-day now that the data exists and the seed-replication
  work is the active thread -- see `cellbin_batch_sigma_slide/`'s section
  above for the full rebuild order.
- `pca_harmony_domain_check.sh`, `domain_check_diagnostics.sh`,
  `domain_check_diagnostics_pending.sh`, `plot_domain_recovery_realmix_prebatch.py`,
  `plot_banksy_domain_panel.py` -- 2026-09-14 domain-check on the real
  (weak-mix) Figure 2/3 data, a sanity check off the strong-mix testbed used
  everywhere else. Writes to `data/figure_3/pca_harmony_domain_figure2_data/`,
  not one of the four dirs above.
- **`plot_domain_vs_celltype_recovery.py`** -- the manuscript panel comparing
  BANKSY domain recovery vs. plain Harmony+Leiden cell-type recovery; stays
  here because it reads across three dirs (`misc/banksy_batch_compare/` +
  `cellbin_batch_sigma_slide/bin/` for domain, `pca_harmony_single_cell/` for
  cell-type), so it doesn't belong to any single one of them.
- `domain_celltype_composition.py`, `plot_composition_recovery.py` -- moved
  from the top level 2026-09-16, superseded: both explicitly called "stale"/
  "superseded" in figure.md (2026-08-27/09-06 entries) once the ARI-based
  domain-vs-celltype framing above was adopted; their outputs are archived at
  `data/figure_3/_archive_20260914/`. Top-level `composition_recovery.py`
  (still essential -- every sweep script's `--h5ad` call feeds the
  `slice_id_leakage_ari` annotation on the manuscript bars) had the other
  half of this same abandoned framing trimmed out 2026-09-16: it used to also
  support a no-args mode scoring all 4 canonical modalities and writing the
  combined `composition_recovery_summary.csv`/`.png` these two scripts read
  -- nothing has called it that way since, so that mode (`score_modality()`,
  `make_plot()`, `--modalities`/`--no-plot`, and `FIG3`/`MODS`) is gone;
  `to_rows()` stays (still used to pretty-print the single result) but
  `--h5ad`/`--out-dir` are now required -- single-file scoring is all this
  script does.
- `run_build_banksy_matrix.sh`, `replot_subsample_umap.py` -- moved from the
  top level 2026-09-16, orphaned: `01_build_banksy_matrix.py` is called
  directly (not via this wrapper) by every current BANKSY script, and
  `replot_subsample_umap.py`'s output dir has never been generated. Neither
  is referenced in figure.md's reproduce commands.
- `plot_banksy_batch_compare.py` -- moved from `misc/banksy_batch_compare/`
  2026-09-16, redundant: a strict subset of `plot_banksy_batch_compare_3mod.py`
  (same cell+spot logic, same JSON paths, just missing bin) -- use the 3-mod
  version, this one has no unique output.

`misc/banksy_batch_compare/run_banksy_batch_compare_umap.sh`/`plot_banksy_batch_compare_umap.py`,
`run_banksy_batch_compare_true_vs_pred.sh`/`plot_banksy_batch_compare_true_vs_pred.py`,
and `cellbin_batch_sigma_slide/run_batch_sigma_slide_umap.sh`/`cell_umap_gap.sh`/
`spot_batch_sigma_slide_weakmix.sh` are real qualitative/validation work (not
duplicates), just not on the critical path to either manuscript panel PNG
(neither reads UMAP diagnostics or the weakmix curve) -- left in place
2026-09-16 by decision, not an oversight.

Full dated narrative for all of the above is in figure.md / the
`figure3_banksy_domain_sweep` memory if you need it -- not duplicated here.

The gene-space Harmony pipeline (`gene_harmony_umap.py` + its `run_*.sh`
wrappers) was deleted 2026-08-27 -- `pca_harmony.py` replaced it and nothing
regenerated its output. `03_clustering_plots.py`'s corresponding
`--pipeline gene_harmony` branch (dead: called into the deleted script) was
stripped 2026-09-16; `--pipeline` stays as a single-choice `pca_harmony` flag
since every caller passes it explicitly. Its own `run_*.sh` step wrappers
(`run_qc_filter.sh`, `run_pca_harmony.sh`, `run_clustering_leiden.sh`,
`run_ari_vs_ground_truth.sh`) moved to `misc/` 2026-09-16 -- see "Steps, if
you need to run them individually" above for why.
