# Count-distribution comparison (Figure 2)

Compares the count distributions of ALBIS simulations with real spatial data
(cells vs Xenium, 8/16 µm bins vs Visium HD, spots vs Visium). The same code
is reused for Figure 5B (ALBIS vs other simulators).

## Layout

```
code/count_distribution/
├── count_distribution.py   all statistics, comparison plots and JSD
├── plot_spatial_spots.py   ALBIS spots drawn at their physical capture radius, one slice
├── _env.sh                 environment (albis-tutorial)
├── smaller_sphere/         Figure 2: one script per simulation-vs-real pairing
└── larger_sphere/          same pairings on a 6 mm sphere (supplementary; not in Figure 2)
```

## What is compared

For two datasets (simulation = `--input`, reference = `--compare-input`),
`count_distribution.py` writes:

| Output | Content |
|---|---|
| `mean_variance_compare.png` | per-gene mean vs variance, Poisson line, NB fit (method-of-moments θ̂) |
| `mean_dropout_compare.png` | per-gene mean vs zero fraction, Poisson and NB predictions |
| `total_counts_compare.png` | counts per observation (density), with JSD |
| `genes_per_cell_compare.png` | detected genes per observation, with JSD |
| `sparsity_summary_compare.png` | % zero entries, % empty observations |
| `raw_norm_log_compare.png` | nonzero values: raw, normalized (10⁴), log1p, each with JSD |
| `comparison_summary.json` | θ̂, medians, sparsity, observation/gene counts |

- **Simulation side:** one slice, stored `slice_id` 4 (the 5th of 10 from the
  bottom; `--all-slices` pools), using post-batch counts (`.X`), like real data.
- **JSD:** Jensen–Shannon divergence, base 2 (0 = identical, 1 = disjoint), on
  60 shared histogram bins; nonzero values are sampled up to 2 × 10⁶ (seed 0).
- **Panel matching** (`--match-panel-size`): the dataset with more genes is
  reduced to the other's gene count by Seurat-v3 highly variable genes. This
  matches panel size, not gene identity.
- `--title-context` changes only the plot titles (default "simulated vs real").

## Figure 2 pairings (`smaller_sphere/`)

| Script | Simulation | Real reference(s) |
|---|---|---|
| `cell_vs_lung_cancer.sh`, `cell_vs_non_diseased_lung.sh` | cell, `log_mu_-2.5_theta_0.40_jitter0.15_bsigma15` | Xenium lung cancer; non-diseased lung |
| `bin8um_same_tissue_as_16um.sh` | bin 8 µm, `packing_pf0p04_bin8um_from16umcfg_log_mu_-2.5_jitter0.6_bsigma07` | Visium HD breast cancer and human pancreas (8 µm) |
| `bin16um_vs_breast_cancer_visium_hd.sh`, `bin16um_vs_human_pancreas_visium_hd.sh` | bin 16 µm, `packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma07` | Visium HD breast cancer; human pancreas (16 µm) |
| `spot_vs_lymph_node_visium.sh`, `spot_vs_tonsil_visium.sh` | spot, `packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_bsigma03` | Visium lymph node; tonsil (probe-based) |

Each script generates the simulation if it is missing, applies QC
(`code/clustering/step00_qc_filter.py`: total > 0, ≥ 3 genes), then runs four modes:

| Mode | Simulation | Panel matching | Use for |
|---|---|---|---|
| `full_panel` | raw | no | mean–variance, mean–dropout (cell) |
| `hvg_matched` | raw | yes | diagnostics |
| `qc_filtered` | QC | no | mean–variance, mean–dropout (bin, spot) |
| `qc_and_hvg_matched` | QC | yes | total counts, genes per observation, sparsity, raw/normalized/log1p |

Bin and spot simulations cover the whole instrument capture window, so their raw
modes are dominated by empty off-tissue observations; use the QC modes for them.

## Reproducing

From `/dcs04/hicks/data/Jan/sim_project`:

```bash
sbatch albis_paper/code/count_distribution/smaller_sphere/<pairing>.sh
```

Outputs: `albis_paper/data/figure_2/smaller_sphere/plots/<pairing>/<mode>/`.
Simulated data: `albis_paper/data/figure_2/smaller_sphere/data/<tag>/`.
Real references (QC'd): `albis_paper/data/real_data_qc/<reference>/`.

Direct use:

```bash
python albis_paper/code/count_distribution/count_distribution.py \
    --modality cell --input <simulation.h5ad> \
    --compare-input <reference.h5ad> --compare-label <name> \
    [--match-panel-size] [--slice-id 4 | --all-slices] --output-dir <dir>
```
