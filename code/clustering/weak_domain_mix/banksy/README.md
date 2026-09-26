# Weak domain mix — BANKSY (batch_sigma 0)

BANKSY on the weak (Figure 2) domain-mix data without batch effect, run at the
λ chosen on the strong mix (`../../strong_domain_mix/banksy/README.md`). λ is
not tuned here: weak-mix domain ARI is ≤ 0.04 at every λ, so a "best" domain λ
would be picking noise.

- Inputs: `data/figure_3/weak_domain_mix/data/<mod>/*_bsigma0_qc.h5ad`
  (`../generate/generate_weakmix_bs0.sh`).
- Why batch_sigma 0: on the batched Figure 2 data, BANKSY at λ > 0 clusters by
  slice (the per-slice batch shift survives neighbour averaging, so Harmony
  has nothing to align). With batch_sigma 0 that collapse disappears.

| File | What |
|---|---|
| `sweep_tasks.tsv` | λ ∈ {0, 0.1, 0.2, 0.3, 0.5} × cell/bin16um/spot (first sweep, 2026-09-24; kept as a record) → `data/figure_3/weak_domain_mix/banksy/sweep/` |
| `run_sweep.sh` | submits rows without results (`DRY_RUN=1` prints only) |
| `domain/`, `cell_type/` | final runs at the strong-mix λ, with UMAP + spatial plots (added once λ is chosen) |

Seeds 101/202 (2026-09-25): data from `../generate/generate_weakmix_seeds.sh` (job 35931529) → `data/figure_3/weak_domain_mix/data/<mod>_seed<n>/`; the usual-σ seeds go next to Figure 2's data (`data/figure_2/smaller_sphere/data/<tag>_seed<n>/`). BANKSY finals 35931538-43, PCA seeds `../pca_harmony/seeds_celltype_panel.sh` (35931544). bin8um stays single-seed (too data-heavy).
