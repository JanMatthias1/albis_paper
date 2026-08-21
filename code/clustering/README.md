# Clustering workflow

Run these from the project root after activating the sim_app tutorial env:

```bash
conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial
python -m pip install -r sim_paper/code/clustering/requirements.txt
python sim_paper/code/clustering/pca_harmony.py --modality spot
python sim_paper/code/clustering/gene_harmony_umap.py --modality spot
python sim_paper/code/clustering/clustering_leiden_louvain.py --modality spot --resolution 0.5
python sim_paper/code/clustering/clustering_leiden_louvain.py --modality spot --algorithm louvain --resolution 0.5
```

Or submit the same steps through SLURM. Each wrapper is a `0-2` array over
`spot`, `bin`, and `cell`; the Leiden/Louvain wrappers run resolutions `0.5`,
`0.1`, and `0.2` within each modality task. All wrappers source `_env.sh` for
the conda-activation boilerplate.

```bash
sbatch sim_paper/code/clustering/run_pca_harmony.sh
sbatch sim_paper/code/clustering/run_gene_harmony_umap.sh
sbatch sim_paper/code/clustering/run_clustering_leiden_louvain.sh
sbatch sim_paper/code/clustering/run_louvain_umap.sh
```

To run `run_pca_harmony.sh` only after two other jobs finish (e.g. to refresh
`pca_harmony` without racing a job that's reading it), use SLURM's own
dependency scheduling instead of adding wait logic to the script:

```bash
sbatch --dependency=afterany:<jobid1>:<jobid2> sim_paper/code/clustering/run_pca_harmony.sh
```

`pca_harmony.py` normalizes/log-transforms the generated `spot`, `bin`, or
`cell` dataset, runs PCA over all genes, runs Harmony on `obs["slice_id"]`,
writes PCA plots before and after Harmony, and plots PC pairs (adjacent pairs
such as PC1 vs PC2, PC3 vs PC4, ...) from the same pre/post-Harmony PCA basis.
The before/after UMAP diagnostic plots run on a subsample (`--umap-max-obs`,
default 50,000; `--no-umap-sample` to disable) since neighbors+UMAP computed
twice at full scale is the slow step for million-cell datasets (e.g. the
`bin` modality's ~6.6M cells) and isn't used downstream — `clustering_leiden_louvain.py`
reads `X_pca_harmony` from the full, non-subsampled output.

`gene_harmony_umap.py` skips PCA for the actual Harmony correction — it runs
Harmony directly on the scaled log-normalized all-gene matrix — but it also
fits a diagnostic PCA independently on the pre- and post-Harmony gene matrices
purely for a PC-pairs sanity check (Harmony was never applied to this basis,
so it's for inspection only, not used downstream). It writes side-by-side UMAP
plots colored by `slice_id`, `domain_true`, and `cell_type_true`, plus the
diagnostic PC-pairs plots.

Both scripts' PC-pairs plots default to coloring by the Harmony `--batch-key`
(`slice_id`); the shared plotting code lives in `pc_pairs.py` (not a
standalone CLI — it's imported by both scripts).

`clustering_leiden_louvain.py` starts from the modality-specific Harmony-corrected PCA file,
uses the first 30 Harmony PCA dimensions, runs Leiden or Louvain clustering
(`--algorithm leiden|louvain`, default `leiden`), computes UMAP, and plots
UMAP colored by `domain_true`, `cell_type_true`, and `cluster_label`. The
output directory/filename tag switches with the algorithm, e.g.
`leiden_res0p5` vs `louvain_res0p5`. `run_louvain_umap.sh` is a thin wrapper
around the same script with `--algorithm louvain` set.

`slice_id` is the right Harmony batch key for the current simulation: sim_app
applies the synthetic gene-expression batch effects per slice and stores the
batch factors in `adata.uns["batch_effect_factors"]`.

## Output layout

Outputs live under `sim_paper/data/clustering/<modality>/<method>/`, where
`<modality>` is `spot`, `bin`, or `cell`, and `<method>` is the step that
produced the output (`pca_harmony`, `gene_harmony_umap`, `leiden_res<tag>`,
`louvain_res<tag>`, and eventually `banksy`). Each method directory holds its `.h5ad` output
alongside a `plots/` subfolder with the PNGs from that step, including a
`plots/pc_pairs/by_<batch_key>/` subfolder for `pca_harmony` and
`gene_harmony_umap` (PC pairs are a diagnostic of those two steps, not a
separate method), e.g.:

```
data/clustering/
  spot/
    pca_harmony/
      simulation_spot_z_pca_harmony.h5ad
      plots/
        pca_before_harmony_by_*.png
        pca_after_harmony_by_*.png
        pc_pairs/
          by_slice_id/
    gene_harmony_umap/
      simulation_spot_z_gene_harmony_umap.h5ad
      plots/
        umap_gene_harmony_before_after_by_*.png
        pc_pairs/
          by_slice_id/
    leiden_res0p5/
      simulation_spot_z_leiden_res0p5.h5ad
      plots/
    louvain_res0p5/
      simulation_spot_z_louvain_res0p5.h5ad
      plots/
  cell/
    pca_harmony/
    gene_harmony_umap/
```

`clustering_leiden_louvain.py` reads its input from `<modality>/pca_harmony/` by default, so
run `pca_harmony.py` for a modality before the other steps.
