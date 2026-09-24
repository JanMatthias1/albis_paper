# Reproducing Figure 3 plots

The active pipeline uses `../manuscript_style.py`: 17 pt titles, 16 pt axis
labels, 12 pt ticks/legend labels, fixed domain/cell-type/modality colors, and
bottom UMAP legends with explicit cluster-to-truth mappings. No script from
`misc/` or `data/figure_3/style_refresh_20260919/` is needed by these workflows.

## Fresh analysis

- `pca_harmony_single_cell/*_celltype_panel.sh`: simulation/QC, plain PCA and
  Harmony, clustering, PCA/UMAP panels, contingency plots, PC-pair diagnostics,
  and per-modality ARI results. Both fresh runs and `--plots-only` use the same
  plotting functions. `replot_typography.sh` can redraw existing outputs.
- `strong_mix/generate_strong_mix_{cell,bin16um,spot}.sh` and
  `cellbin_batch_sigma_slide/batch_sigma_slide_seeds.sh`: simulation/QC,
  BANKSY/Harmony, resolution selection and scores, followed automatically by
  `plot_banksy_results.py` for PCA, before/after UMAP, matched-cluster UMAP and
  contingency plots. Plotting adds display-embedding runtime to each job.
- `spatial_deconvolution/{weak_domain,strong_domain}/run_rctd_spot*.sh`: fitting followed automatically by
  `plot_rctd_results.py` for spatial zooms and the error panel. Independent-seed
  fitting requires the separately generated reference datasets, as before.
- After the analysis jobs finish, run `bash plot_figure3_summaries.sh` for the
  cross-modality ARI bars and batch-strength curve. Summary scripts retain their
  existing missing-result behavior; a complete curve requires all intended jobs.
- For the averaged RCTD results, run `spatial_deconvolution/average_rctd_seeds.py
  --config weak_mix` (and `--config strong_mix`), then
  `spatial_deconvolution/plot_rctd_results.py` with the `weak_mix_seed_avg`
  (`strong_mix_seed_avg`) directory as its positional argument.

Submit the analysis shell scripts through Slurm; they can be expensive.
The summary launcher is a local plotting-only command.

## Replot BANKSY without rerunning analysis

Using the `sim_paper/env/albis-tutorial/bin/python` interpreter:

```bash
python sim_paper/code/clustering/plot_banksy_results.py \
  --input /path/to/banksy_matrix/simulation_spot_z_banksy_pca_harmony_qc.h5ad \
  --cluster-input /path/to/ari/simulation_spot_z_ari_recovery.h5ad
```

Writes into each input's sibling `plots/` directory. Supply `--output-dir` to
render into an alternate root. Only observation metadata and embeddings are
read; analysis h5ads are never modified. The default display subsample is 50,000
observations, seed 0, 15 neighbors, matching the original refresh. UMAP caches
are validated against observation order and PCA values before reuse. On a fresh
run they are computed directly; old style-refresh caches are not required.

Styling is reproducible through the active code. Pixel-identical fresh
simulation/embedding results also depend on the same input data, seeds and
software environment; this update does not rerun or revalidate the full
scientific workflow. Historical/diagnostic scripts under `misc/` are excluded.
