# Figure 2 methodology: sim-vs-real count-distribution comparison

Snapshot of the current approach for the Figure 2 panels (statistical/
count-distribution fidelity of simulated data vs. real reference data), so
this can be picked up later for fine-tuning without re-deriving the
reasoning. Parameter values here are current best guesses, not final --
expect them to change. For the full tuning history (sweeps tried, numbers
at each step, open leads), see `figure.md` and `DATA_VERSIONS.md` at the
`sim_paper/` root; this doc only covers the comparison *methodology*, not
the calibration history.

## What's being compared

Seven modality-vs-real pairings, one script each under this directory
(`<pairing>.sh` and `<pairing>_with_batch.sh`):

| Modality | Real reference | Sim config tag |
|---|---|---|
| cell | Xenium `non_diseased_lung` | `log_mu_-2.5_theta_0.25_jitter0.15` |
| cell | Xenium `lung_cancer` | `log_mu_-2.5_theta_0.25_jitter0.15` (same) |
| bin (8um) | Visium HD `breast_cancer_visium_hd` | `packing_pf0p04_log_mu_0.0` |
| bin (8um) | Visium HD `human_pancreas_visium_hd` | `packing_pf0p04_log_mu_0.0` (same) |
| bin (16um) | Visium HD `breast_cancer_visium_hd_16um` | `packing_pf0p04_bin16um_log_mu_-2.5` |
| bin (16um) | Visium HD `human_pancreas_visium_hd_16um` | `packing_pf0p04_bin16um_log_mu_-2.5` (same) |
| spot | Visium `breast_cancer_visium` | `packing_pf0p04_log_mu_-2.5` |

Each script is self-contained: generates the sim data if missing, QC-filters
it (`code/clustering/00_qc_filter.py` -- drops fully-empty and near-empty
observations, a graph-fragmentation fix borrowed from the Figure 3
clustering pipeline, not a real-data-parity correction), then runs
`count_distribution.py` in four modes (see below). Bin/16um real references
are resampled from the same raw 10x `binned_outputs/` at 16um instead of
8um (`code/real_data_qc/misc/run_visium_hd_qc_16um.sh`) -- both resolutions
are canonical Visium HD configurations, not a replacement of one by the
other.

## Two methodological choices, currently run in parallel

**1. Single representative slice, not pooled.** Every sim dataset has 10
z-plane slices through the simulated tissue sphere (`slice_id` 0-9). Real
reference data is a single tissue section. Pooling all 10 slices is not
apples-to-apples, and empirically two of them (`slice_id` 0 and 9, nearest
the sphere's poles) are geometrically degenerate -- their cross-sectional
disc is much smaller than at the equator, so most observations land
off-tissue (~21% empty bins vs. <0.5% for the 8 interior slices on the
bin16um config). Pooling silently dragged the aggregate stats down.
**Decision: restrict to `slice_id=5`** (near-equator, representative)
via `count_distribution.py --slice-id 5`, applied to the sim/primary side
only (real data has no `slice_id`).

**2. Batch effect: kept both ways, not yet decided.** The sim's synthetic
per-slice technical batch shift (`batch_sigma`, used elsewhere to test
Harmony correction) can be included or excluded from the count-distribution
stats via `count_distribution.py --use-batch-effect`. Default excludes it
(uses `counts_pre_batch`); the flag makes it use `X` uniformly instead.
Real data has no "pre-batch" version -- whatever technical noise a real
section carries is just baked into its counts -- so arguably including the
sim's batch effect is the fairer comparison now that we're down to one
slice. It is *not* a small effect (bin16um `theta_hat` at slice 5: 1.22
pre-batch vs. 0.18 post-batch), so this is a real methodological choice,
not a rounding correction.

Both variants are generated side by side so neither is lost:
- `sim_paper/data/count_distribution/figure_2/` -- slice 5, pre-batch
  (`counts_pre_batch`, the historical default).
- `sim_paper/data/count_distribution/figure_2_with_batch/` -- slice 5,
  post-batch (`X`, includes the synthetic batch shift).

**Not yet decided which becomes the one cited in the manuscript.** Compare
both once landed; `figure_2_with_batch` is the more defensible
apples-to-apples choice on paper, but check the actual plots before
committing.

## The four comparison modes

Each pairing script runs `count_distribution.py` four times, into
`<OUT_ROOT>/{full_panel,hvg_matched,qc_filtered,qc_and_hvg_matched}/`:

| Mode | Sim input | `--match-panel-size` | Use for |
|---|---|---|---|
| `full_panel` | raw (not QC'd) | no | `mean_variance`/`mean_dropout` dispersion fit -- full gene panel, unfiltered; HVG-selection or QC would bias the dispersion estimate itself. |
| `hvg_matched` | raw (not QC'd) | yes | `sparsity_summary`/empty-row-rate -- panel-matched but not QC'd, since QC-filtering sim here would hide the zero-inflation these stats exist to report. |
| `qc_filtered` | QC'd | no | `total_counts` (secondary use) with QC'd sim, full gene panel. |
| `qc_and_hvg_matched` | QC'd | yes | `total_counts`/`genes_per_cell` -- both corrections together, the fair like-for-like for these two panels specifically. |

`--match-panel-size` HVG-subsets whichever side (sim or real) has the
larger gene panel down to the smaller side's count, using `seurat_v3`
highly-variable-gene selection (not a random subsample) so the comparison
isn't just measuring panel-size differences. It's symmetric: usually
shrinks the real reference (Visium(HD) has 18k-36k genes vs. sim's 556),
but for cell vs. Xenium sim's 556-gene panel is actually the larger one
(Xenium has only 392), so it shrinks sim instead.

Only cite the mode each panel type is designed for -- e.g. never read
`mean_variance_compare.png` from `qc_and_hvg_matched`, always from
`full_panel`.

## Rerunning

```bash
cd /dcs04/hicks/data/Jan/sim_project
sbatch sim_paper/code/count_distribution/<pairing>.sh              # pre-batch
sbatch sim_paper/code/count_distribution/<pairing>_with_batch.sh   # post-batch
```

Each is idempotent for the generate/QC steps (skips if the sim `.h5ad`
already exists) but always reruns all 4 `count_distribution.py` modes.
To change the sim config, edit `SIM_TAG` + the `generate_simulation_noisy.py`
flags together at the top of the relevant script(s) and resubmit -- doesn't
touch the other pairings.

To point a run at a different slice or toggle the batch effect manually
(e.g. for a quick one-off check), the two relevant `count_distribution.py`
flags are `--slice-id <int>` and `--use-batch-effect` (see
`count_distribution.py --help` or `README.md` for the rest of the CLI).

## Known open gaps (not addressed by this methodology, tracked in figure.md)

- Pancreas is a known-bad calibration target for bin (extreme real
  overdispersion from a few dominant hormone genes) -- not expected to
  match well regardless of slice/batch choice.
- Spot's `genes_per_cell` still mismatches real even panel-matched --
  open issue, unrelated to slice/batch.
- QC-parity gap: real bin/spot data goes through SpotSweeper local-outlier
  QC before comparison; sim only gets the minimal `00_qc_filter.py` pass
  (empty/near-empty only). Parked, not addressed here.
- Per-modality `batch_sigma` values (used for Figure 3's Harmony-correction
  demo) are a separate, larger tuning question from the include/exclude
  choice in `--use-batch-effect` here -- this doc doesn't change what
  `batch_sigma` each sim config used to *generate* the data, only whether
  Figure 2's stats read the pre- or post-batch layer of what's already
  there.
