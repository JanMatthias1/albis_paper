# Count distribution diagnostics

Run from the project root after activating the sim_app tutorial env (same env
as `sim_paper/code/clustering`, no extra packages needed):

```bash
conda activate /dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial
python sim_paper/code/count_distribution/count_distribution.py --modality spot
python sim_paper/code/count_distribution/count_distribution.py --modality bin
python sim_paper/code/count_distribution/count_distribution.py --modality cell
```

Or submit all three modalities through SLURM (a `0-2` array over `spot`,
`bin`, `cell`, sourcing `_env.sh` for the conda-activation boilerplate):

```bash
sbatch sim_paper/code/count_distribution/run_count_distribution.sh
```

`count_distribution.py` produces four PNGs per modality under
`sim_paper/data/count_distribution/<modality>/`:

- `mean_variance.png` -- per-gene mean vs variance (log-log), with the
  Poisson line (`var = mean`) and a single method-of-moments NB line
  (`var = mean + mean^2/theta`) overlaid.
- `mean_dropout.png` -- per-gene mean vs observed zero fraction, with the
  Poisson- and NB-predicted zero-probability curves overlaid (same
  `theta` as `mean_variance.png`). Tests whether the NB fit from the first
  two moments also explains the zero rate, or whether the data needs
  zero-inflation on top of NB.
- `total_counts.png` -- per-cell library-size (total counts) histogram.
- `raw_norm_log.png` -- histogram of nonzero matrix entries at three
  pipeline stages: raw counts, target-sum normalized, and
  `log1p(normalized)`. Uses the same `normalize_total`/`log1p` steps as
  `clustering.py`, on a random sample of up to `--sample-size` nonzero
  entries.

`mean_variance.png`, `mean_dropout.png`, and `total_counts.png` are
computed from `adata.layers["counts_pre_batch"]` when present (the NB draw
before the synthetic batch-effect multiplier), falling back to `adata.X`
otherwise -- this checks the underlying generative count model, not the
batch-perturbed data. `raw_norm_log.png` always uses `adata.X`, since
that's what the clustering pipeline actually consumes.

## Comparing against real datasets

`--input`, `--output-dir`, and `--modality` are independent: point `--input`
at any `.h5ad` with raw counts in `.X`, pick an `--output-dir` (or reuse
`--modality` as a label to land under `sim_paper/data/count_distribution/`),
and the same three plots are produced. `--modality` is not restricted to
`spot`/`bin`/`cell` -- it's only used to build default paths and plot
titles. Real data without a `counts_pre_batch` layer just falls back to
`adata.X` for the mean-variance/total-counts panel.

```bash
python sim_paper/code/count_distribution/count_distribution.py \
    --input /path/to/real_dataset.h5ad \
    --modality xenium_breast
```

## Overlaying sim vs. real on the same plots

Pass `--compare-input` (plus `--compare-label`) to overlay a second dataset on
the same four diagnostics instead of producing separate PNGs per dataset --
useful for checking whether the simulated count distribution actually looks
like real Xenium data:

```bash
python sim_paper/code/count_distribution/count_distribution.py \
    --modality cell \
    --compare-input sim_paper/data/real_data_qc/non_diseased_lung/non_diseased_lung_qc.h5ad \
    --compare-label non_diseased_lung
```

This produces, under `sim_paper/data/count_distribution/<modality>_vs_<compare-label>/`:

- `mean_variance_compare.png`, `mean_dropout_compare.png` -- both datasets'
  per-gene scatter clouds and NB fits overlaid on one plot (one shared
  Poisson reference line).
- `total_counts_compare.png` -- density-normalized (not raw-count)
  histograms of total counts per cell, so datasets with very different
  numbers of cells are still comparable by shape.
- `raw_norm_log_compare.png` -- same density-normalized overlay, for each of
  the raw/normalized/log1p stages.

Colors are fixed by dataset role (`--modality` = blue, `--compare-label` =
orange) across all four plots, not by draw order.
