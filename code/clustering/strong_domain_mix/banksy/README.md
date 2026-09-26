# Strong domain mix — BANKSY λ (batch_sigma 0)

Chooses BANKSY's λ per modality and target (domain, cell type) on the strong
domain-mix data. The weak mix is then run at the **same** λ
(`../../weak_domain_mix/banksy/`), so the two rows of
`banksy_vs_pca_recovery.png` differ only in the data.

- Inputs: `data/figure_3/strong_domain_mix/batch_sigma_slide/<mod>/bs0[.0]/*_qc.h5ad`
  (the batch-σ sweep's no-batch point; nothing regenerated here).
- k_geom: cell 60, bin16um 100, spot 8 (same as the batch-σ sweep).
- Each run scores both `domain_true` and `cell_type_true`
  (`../../run_banksy_lambda.sh`), one seed.

| File | What |
|---|---|
| `sweep_tasks.tsv` | λ ∈ {0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0} × cell/bin16um/spot → `data/figure_3/strong_domain_mix/banksy/sweep/<mod>/lam<λ>/` |
| `run_sweep.sh` | submits rows without results (`DRY_RUN=1` prints only) |
| `domain/`, `cell_type/` | final runs at the chosen λ, with UMAP + spatial plots (added once λ is chosen) |

## Status

- 2026-09-24: λ 0–0.5 run (first sweep, then in `data/banksy_cell/`).
- 2026-09-25: moved here; λ 0.8 and 1.0 added because cell and bin16um domain
  ARI was still rising at 0.5 (jobs 35926530 cell, 35926531 bin16um, 35926532 spot).
  bin16um at λ 0.8 and 1.0 clusters by slice (clusters vs slice_id ARI 0.999,
  Leiden stuck at 10 clusters down to resolution 0.01) even though this data
  has no batch effect; at λ 0.5 it is 0.015 vs slice and 0.80 vs domain.
  Cause not investigated.

## Sweep results (ARI; * = Leiden missed the target cluster count)

| | λ 0 | 0.1 | 0.2 | 0.3 | 0.5 | 0.8 | 1.0 |
|---|---|---|---|---|---|---|---|
| cell domain | 0.05 | 0.04 | 0.11 | 0.40* | 0.44 | 0.44 | 0.55 |
| bin16um domain | 0.00 | 0.17 | 0.54 | 0.77 | 0.80 | 0.02* | 0.02* |
| spot domain | 0.33 | 0.48 | 0.60 | 0.77 | 0.58 | 0.78 | 0.68 |
| cell cell type | 0.65 | 0.66 | 0.67 | 0.37 | 0.07 | 0.06 | 0.06 |
| bin16um cell type | 0.00 | 0.24 | 0.11 | 0.11 | 0.10 | 0.00* | 0.00* |
| spot cell type | 0.35 | 0.34 | 0.32 | 0.30 | 0.29 | 0.28 | 0.27 |

## Chosen λ

Decided 2026-09-25 (used for both mixes; `domain/final_tasks.tsv`,
`cell_type/final_tasks.tsv`):

| | cell | bin16um | spot |
|---|---|---|---|
| domain | 1.0 | 0.5 | 0.3 |
| cell type | 0.2 | 0.1 | 0 |

- cell domain λ 1.0 is the top of the range (neighbourhood only). Accepted for
  domains: a single cell carries one cell type, so its domain (a cell-type
  mix) is only visible in its neighbourhood. The λ 0.5 / 0.8 / 1.0 gap
  (0.44 / 0.44 / 0.55) is near run-to-run noise: the batch-σ plot's identical
  λ 0.5 run scored 0.48 (embeddings differ by ≤ 0.007, Leiden flips), and the
  three simulation seeds at batch_sigma 0 span 0.48–0.62.
- spot domain: λ 0.3 (0.77) and 0.8 (0.78) tie; 0.3 chosen as the smaller λ.
- bin16um domain stays 0.5: at λ ≥ 0.8 it clusters by slice (see Status).

Final runs rerun BANKSY (same flags) with UMAP + spatial plots, so their ARI
can differ from the sweep by that run-to-run noise.

Run: `bash domain/run_final.sh`, `bash cell_type/run_final.sh` (and the same in
`../../weak_domain_mix/banksy/`) → `data/figure_3/<mix>/banksy/{domain,cell_type}/<mod>/lam<λ>/`.

The batch-σ plot (`../batch_sigma_slide/`) keeps its original domain λ
(cell 0.5, bin16um 0.5, spot 0.1).

## Simulation seeds 101/202 (submitted 2026-09-25)

Every Figure 3 bar is now the mean ± SD over simulation seeds 2025/101/202
(`final_tasks.tsv` has a `seed` column; seed folders end in `_seed<n>`).
Jobs: weak-mix generation 35931529 (`../../weak_domain_mix/generate/generate_weakmix_seeds.sh`);
strong BANKSY finals 35931530-35 (inputs: batch_sigma_slide/<mod>/bs0_seed<n>/);
strong PCA seeds 35931536 (strongmix_celltype_panel.sh --array=6-17);
weak BANKSY finals 35931538-43 and weak PCA seeds 35931544 (afterany 35931529).
bin16um slice-split diagnostic: sweep row 21 (λ 0.8, seed 101, KEEP_H5AD=1), job 35931537.
