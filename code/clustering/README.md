# Clustering and deconvolution (Figure 3)

Code for spatial-domain and cell-type recovery on ALBIS simulations
(cells, 16 µm bins, spots), the batch-effect sweep, and RCTD deconvolution.
Outputs are written to `sim_paper/data/figure_3/`, which mirrors this layout.

## Layout

```
code/clustering/
├── step00_qc_filter.py               drop empty / near-empty observations (total > 0, ≥ 3 genes)
├── step01_pca_harmony.py             genes → normalize → log1p → PCA (→ Harmony)
├── step01_build_banksy_matrix.py     BANKSY spatial embedding → PCA (→ Harmony)
├── step02_leiden_resolution_sweep.py Leiden, resolution searched to the true group count; ARI
├── step03_cluster_and_plot.py        UMAP / contingency plots at that resolution
├── composition_recovery.py           ARI and slice-leakage scores for the batch-σ sweep
├── pc_pairs.py, plot_banksy_results.py   diagnostic plots
├── run_banksy_lambda.sh, submit_banksy_lambda.sh   BANKSY at a fixed λ from a task table
├── run_figure3_global.sh             submits the Figure 3 stages with dependencies
├── plot_figure3_summaries.sh         redraws the summary figures
├── no_harmony_domain/                batch σ = 0 domain recovery, no Harmony (λ sweep → final runs)
├── no_harmony_cell_type/             batch σ = 0 cell-type recovery, no Harmony (λ sweep → final runs)
├── strong_domain_mix/
│   ├── generate/                     strong-mix simulations (+ check against Figure 2 settings)
│   ├── batch_sigma_slide/            batch-σ sweep: BANKSY + Harmony, 3 seeds, final plot
│   ├── pca_harmony/                  cell-type panels on strong-mix data
│   └── banksy_harmony_batch_zero/    BANKSY + Harmony at σ = 0 (reference only)
├── weak_domain_mix/
│   ├── generate/                     σ = 0 and extra-seed versions of the Figure 2 data
│   ├── pca_harmony/                  cell-type panels on the Figure 2 data
│   └── banksy_harmony_batch_zero/    BANKSY + Harmony at σ = 0 (reference only)
├── ari_recovery_summary/             recovery bar figures
└── spatial_deconvolution/            RCTD on spots (weak_domain/, strong_domain/), plots
```

## Pipeline

```
QC (step00) → embedding (step01: genes or BANKSY) → Leiden at the true group count (step02) → ARI
```

Harmony is applied only when a batch effect is simulated:

| Setting | Embedding | Harmony |
|---|---|---|
| batch σ = 0 (`no_harmony_*`) | BANKSY or genes → PCA → Leiden | no |
| tuned batch σ (cell 1.5, bin16um 0.7, spot 0.3) | genes → PCA → Harmony → Leiden | yes |
| batch-σ sweep (`batch_sigma_slide/`) | BANKSY → PCA → Harmony → Leiden | yes |

At σ = 0 the batch factor is 1, but counts pass through the same Poisson
resampling as batched data, so technical noise is matched across conditions.

## Parameters

| | cell | bin16um | spot |
|---|---:|---:|---:|
| BANKSY k_geom | 60 | 100 | 8 |
| BANKSY λ, domain (σ = 0 and batch-σ sweep) | 0.3 | 0.3 | 0.5 |
| Tuned batch σ | 1.5 | 0.7 | 0.3 |

Batch-zero λ were selected on the strong mix, seed 2025 (highest ARI among runs
with the true cluster count), and reused for the weak mix and seeds 101, 202;
see `data/figure_3/no_harmony_*/selected_parameters.json`. Seed 2025 is also
included in the reported means.

## Reproducing

From `/dcs04/hicks/data/Jan/sim_project`:

```bash
# cell-type panels, RCTD, batch-σ sweep, summary plots
bash sim_paper/code/clustering/run_figure3_global.sh

# batch-zero, no-Harmony experiments (each: sweep → select → final → summary)
sbatch --array=0-20 sim_paper/code/clustering/no_harmony_domain/run.sh sweep
sbatch sim_paper/code/clustering/no_harmony_domain/run.sh select
sbatch --array=0-35 sim_paper/code/clustering/no_harmony_domain/run.sh final
sbatch sim_paper/code/clustering/no_harmony_domain/run.sh summary
# same four stages for no_harmony_cell_type/run.sh (after no_harmony_domain)

# figures only
bash sim_paper/code/clustering/plot_figure3_summaries.sh
```

Environments (`sim_paper/env/`): `albis-tutorial` (`_env.sh`; QC, PCA/Harmony,
Leiden, plots), `albis-banksy` (`_env_banksy.sh`; BANKSY), `rctd` (RCTD, R).
Python requirements: `requirements.txt`.

Single-slice plots use stored `slice_id` 4 (the 5th of 10 slices from the bottom).
