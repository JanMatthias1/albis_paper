# Clustering workflow

Run these from the project root after activating the sim_app tutorial env:

```bash
conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial
python -m pip install -r sim_paper/code/clustering/requirements.txt
python sim_paper/code/clustering/clustering.py
python sim_paper/code/clustering/plot_pc_pairs.py --n-pcs 12
python sim_paper/code/clustering/leiden_umap.py --resolution 0.5
```

Or submit the same steps through SLURM:

```bash
sbatch sim_paper/code/clustering/run_clustering.sh
sbatch sim_paper/code/clustering/run_leiden_umap.sh --resolution 0.5
```

`clustering.py` normalizes/log-transforms the generated Visium spot dataset, runs
PCA over all genes, runs Harmony on `obs["slice_id"]`, and writes PCA plots
before and after Harmony.

`plot_pc_pairs.py` reads the PCA/Harmony file and writes pre-Harmony and
post-Harmony pairwise PC scatter plots, such as PC1 vs PC2, PC1 vs PC3, and so
on. It defaults to all PCs and `slice_id` coloring; use `--n-pcs`, `--color-key`,
or `--no-sample` to tune the output.

`leiden_umap.py` starts from the Harmony-corrected PCA file, uses the first 30
Harmony PCA dimensions, runs Leiden clustering, computes UMAP, and plots UMAP
colored by `domain_true`, `cell_type_true`, and `cluster_label`.

`slice_id` is the right Harmony batch key for the current simulation: sim_app
applies the synthetic gene-expression batch effects per slice and stores the
batch factors in `adata.uns["batch_effect_factors"]`.
