# Clustering workflow

Run these from the project root after activating the sim_app tutorial env:

```bash
conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial
python -m pip install -r sim_paper/code/clustering/requirements.txt
python sim_paper/code/clustering/clustering.py --modality spot
python sim_paper/code/clustering/plot_pc_pairs.py --modality spot --n-pcs 12
python sim_paper/code/clustering/gene_harmony_umap.py --modality spot
python sim_paper/code/clustering/leiden_umap.py --modality spot --resolution 0.5
```

Or submit the same steps through SLURM. Each wrapper is a `0-2` array over
`spot`, `bin`, and `cell`; the Leiden wrapper runs resolutions `0.5`, `0.1`,
and `0.2` within each modality task.

```bash
sbatch sim_paper/code/clustering/run_clustering.sh
sbatch sim_paper/code/clustering/run_plot_pc_pairs.sh
sbatch sim_paper/code/clustering/run_gene_harmony_umap.sh
sbatch sim_paper/code/clustering/run_leiden_umap.sh
```

`clustering.py` normalizes/log-transforms the generated `spot`, `bin`, or
`cell` dataset, runs PCA over all genes, runs Harmony on `obs["slice_id"]`, and
writes PCA plots before and after Harmony.

`plot_pc_pairs.py` reads the PCA/Harmony file and writes two grid figures: one
pre-Harmony and one post-Harmony. Each figure shows adjacent PC pairs, such as
PC1 vs PC2, PC3 vs PC4, and so on. It defaults to all PCs and `slice_id`
coloring; use `--n-pcs`, `--color-key`, or `--no-sample` to tune the output.

`gene_harmony_umap.py` skips PCA, runs Harmony directly on the scaled
log-normalized all-gene matrix, computes UMAP before and after Harmony, and
writes side-by-side UMAP plots colored by `slice_id`, `domain_true`, and
`cell_type_true`.

`leiden_umap.py` starts from the modality-specific Harmony-corrected PCA file,
uses the first 30 Harmony PCA dimensions, runs Leiden clustering, computes
UMAP, and plots UMAP colored by `domain_true`, `cell_type_true`, and
`cluster_label`.

`slice_id` is the right Harmony batch key for the current simulation: sim_app
applies the synthetic gene-expression batch effects per slice and stores the
batch factors in `adata.uns["batch_effect_factors"]`.
