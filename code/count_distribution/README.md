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

`count_distribution.py` produces three PNGs per modality under
`sim_paper/data/count_distribution/<modality>/`:

- `mean_variance.png` -- per-gene mean vs variance (log-log), with the
  Poisson line (`var = mean`) and a single method-of-moments NB line
  (`var = mean + mean^2/theta`) overlaid.
- `total_counts.png` -- per-cell library-size (total counts) histogram.
- `raw_norm_log.png` -- histogram of nonzero matrix entries at three
  pipeline stages: raw counts, target-sum normalized, and
  `log1p(normalized)`. Uses the same `normalize_total`/`log1p` steps as
  `clustering.py`, on a random sample of up to `--sample-size` nonzero
  entries.

`mean_variance.png` and `total_counts.png` are computed from
`adata.layers["counts_pre_batch"]` when present (the NB draw before the
synthetic batch-effect multiplier), falling back to `adata.X` otherwise --
this checks the underlying generative count model, not the batch-perturbed
data. `raw_norm_log.png` always uses `adata.X`, since that's what the
clustering pipeline actually consumes.

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
