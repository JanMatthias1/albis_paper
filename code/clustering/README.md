# Clustering workflow (Figure 3)

Two conda envs are involved:

```bash
# Everything except the BANKSY build step
conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial
python -m pip install -r sim_paper/code/clustering/requirements.txt

# BANKSY build step only (01_build_banksy_matrix.py) -- banksy_py pins an
# older scanpy/numpy/anndata/pandas/scikit-learn/scipy stack, incompatible
# with sim-app-tutorial. See sim_app/env/create_banksy_env.sh.
conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-banksy
```

All scripts here source `_env.sh` (sim-app-tutorial) or `_env_banksy.sh`
(sim-app-banksy) for the conda-activation boilerplate.

## Reproduce Figure 3: run one panel script per modality x task

The fastest path to reproducing a Figure 3 panel is one of these six
self-contained scripts, one per (modality x cell-type-or-domain panel):

```bash
sbatch sim_paper/code/clustering/spot_celltype_panel.sh
sbatch sim_paper/code/clustering/bin_celltype_panel.sh
sbatch sim_paper/code/clustering/cell_celltype_panel.sh
sbatch sim_paper/code/clustering/spot_domain_panel.sh
sbatch sim_paper/code/clustering/bin_domain_panel.sh
sbatch sim_paper/code/clustering/cell_domain_panel.sh
```

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
`--plots-only` flag (`pca_harmony.py`, `clustering_leiden_louvain.py`)
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
- **Domain panels** (`*_domain_panel.sh`): BANKSY -> PCA -> Harmony ->
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
     key sim_app's synthetic batch effects are applied per). Writes PCA and
     before/after-Harmony UMAP diagnostic plots plus PC-pairs plots
     (`pc_pairs.py`, imported, not a standalone CLI).
     ```bash
     python sim_paper/code/clustering/pca_harmony.py --modality bin --input <qc.h5ad> --output <out.h5ad>
     ```
   - `01_build_banksy_matrix.py` (needs the `sim-app-banksy` env): BANKSY ->
     PCA -> Harmony, staggering each slice's coordinates first (all slices
     share one local x/y range, so an unstaggered neighbor graph would treat
     different slices as spatial neighbors). Writes `obsm["X_pca_harmony"]`
     in the same shape `pca_harmony.py` does, so `clustering_leiden_louvain.py`
     needs zero changes to consume either.
     ```bash
     python sim_paper/code/clustering/01_build_banksy_matrix.py --modality cell --packing-tag <tag>
     ```
4. **`clustering_leiden_louvain.py`** -- Leiden or Louvain clustering
   (`--algorithm`, default `leiden`) + UMAP from either pipeline's output
   (`--pipeline pca_harmony|gene_harmony`, default `pca_harmony`), plots UMAP
   colored by `domain_true`/`cell_type_true`/`cluster_label`, and (if both a
   `--cluster-key` and a `*_true` ground-truth column are present) a
   ground-truth-vs-predicted UMAP + contingency heatmap.
   ```bash
   python sim_paper/code/clustering/clustering_leiden_louvain.py --modality bin --resolution 0.5
   ```
5. **`ari_vs_ground_truth.py`** -- Figure 3B: binary-searches the Leiden
   resolution until the cluster count matches the true category count (so
   ARI isn't confounded by over/under-clustering), separately for
   `domain_true` and `cell_type_true`.
   ```bash
   python sim_paper/code/clustering/ari_vs_ground_truth.py --modality bin --packing-tag <tag>
   ```
6. **`plot_ari_recovery.py`** -- grouped bar charts of ARI recovery across
   modality x pipeline, auto-discovering every
   `data/clustering_<tag>/*/*ari_recovery*/ari_summary_<modality>.json` it
   can find.
   ```bash
   python sim_paper/code/clustering/plot_ari_recovery.py --packing-tag <tag>
   ```

Corresponding `run_*.sh` SLURM wrappers exist for steps 2-6 (`0-2` array
over spot/bin/cell where relevant) for running a step standalone/chained via
`--dependency=afterok:<jobid>` instead of through a panel script -- see
figure.md's "Reproduce" code blocks for worked examples.

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

**Domain panels are not yet under `data/figure_3/`** -- known inconsistency,
not yet resolved (flagged 2026-08-26). Each modality's domain-panel dataset
carries its own tag (differs because `batch_sigma`/`domain_type_mix` differ
per modality) and lives at `data/clustering_<tag>/<modality>/`, mirroring
the same three-step layout (`banksy_pca_harmony_qc/`, `banksy_ari_recovery/`,
`leiden_pca_banksy_domain_matched/`) under a different root:

| modality | domain-panel tag |
|---|---|
| bin | `packing_pf0p04_bsigma05_strongdomainmix` |
| cell | `log_mu_-2.5_theta_0.25_strongdomainmix` |
| spot | `packing_pf0p04_bsigma015_strongdomainmix` |

Check the modality's `*_domain_panel.sh` header for its current `SIM_TAG` --
these get retuned as work continues, so treat the table above as a snapshot,
not a guarantee.

## `misc/`

Two different kinds of thing live here, not one:

- **Superseded/deprecated**: `gene_harmony_umap.py`, `run_gene_harmony_umap.sh`,
  `run_louvain_umap.sh` (the gene-space Harmony pipeline `pca_harmony.py`
  replaced; `clustering_leiden_louvain.py --pipeline gene_harmony` still
  supports reading its output if regenerated, but nothing currently
  regenerates it), `pc_pairs.py` used to be here too -- **moved back** to
  top-level 2026-08-26 after discovering `pca_harmony.py` still imports it
  directly (it's a live dependency, not dead code; same for the `run_*.sh`
  step wrappers below, which figure.md's own documented reproduce commands
  reference at their top-level path -- don't re-move things into `misc/`
  without grepping `figure.md`/`DATA_VERSIONS.md` for existing references
  first).
- **One-off exploration/tuning scripts**, kept for provenance:
  `tune_banksy_bin_kgeom.sh`, `tune_banksy_cell.sh` (BANKSY hyperparameter
  sweeps, see figure.md), `replot_*.sh` (one-off re-plots after a
  plot-formatting fix, superseded once the fix is merged into the scripts
  themselves).
