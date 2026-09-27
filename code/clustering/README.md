# Clustering workflow (Figure 3)

Current workflow, updated 2026-09-18. Cell-type recovery and deconvolution
reuse Figure 2 smaller-sphere data. Domain recovery and batch-effect sweeps
use separate simulations generated with `--strong-domain-mix`.

## Layout (2026-09-25)

Scripts are grouped by dataset; `data/figure_3/` mirrors the same tree.

```
code/clustering/
├── 00_qc_filter.py … step03_cluster_and_plot.py, plot_banksy_results.py   shared pipeline steps
├── run_banksy_lambda.sh, submit_banksy_lambda.sh                         BANKSY at a fixed λ (task-table driven)
├── strong_domain_mix/
│   ├── generate/           generate_strong_mix_{cell,bin16um,spot}.sh, check_matches_figure2.py
│   ├── batch_sigma_slide/  batch-σ sweep: seeds, final plot, cliff UMAPs
│   ├── banksy/             λ sweep at batch_sigma 0 (+ domain/, cell_type/ final runs)
│   └── pca_harmony/        strongmix_celltype_panel.sh (bs0 + usual σ)
├── weak_domain_mix/
│   ├── generate/           generate_weakmix_bs0.sh (batch_sigma 0 only; batched = Figure 2 data)
│   ├── banksy/             λ sweep at batch_sigma 0 (+ domain/, cell_type/ final runs)
│   └── pca_harmony/        {cell,bin,bin16um,spot}_celltype_panel.sh (Figure 2 data), bs0_celltype_panel.sh
├── ari_recovery_summary/   plots across mixes/pipelines
└── spatial_deconvolution/  RCTD (weak_domain/, strong_domain/)
```

Moved from `strong_mix/`, `cellbin_batch_sigma_slide/`, `pca_harmony_single_cell/`
and `data/banksy_cell/`; see `data/figure_3/README.md` for the data side. The old
`cellbin_batch_sigma_slide` path is retired; Figure 4B STAGATE and RCTD read the new one.

## Pipeline overview

`02_leiden_resolution_sweep.py` and `step03_cluster_and_plot.py` are
**siblings**, not a chain -- both are called with the identical embedding
`--input` (see any `*_celltype_panel.sh`); `step03` does not read `02`'s
output h5ad, it only reuses the resolution *value* `02`'s binary search
already found. `composition_recovery.py` is the domain-sweep branch's
counterpart to `step03`'s plotting branch -- both sit at the same depth,
consuming `02`'s resolution-matched output, just for different purposes
(qualitative UMAP/heatmap panels vs. numeric ARI+leakage scoring). This is
why `composition_recovery.py` is not named `step04_...`: it isn't downstream
of `step03`, it's a fork alongside it.

```
00_qc_filter.py
     |
     +-- 01.2_pca_harmony.py       (plain)   --+
     +-- 01_build_banksy_matrix.py (BANKSY)  --+--> obsm["X_pca_harmony"]
                                                          |
                    +-------------------------------------+-------------------------------------+
                    |                                                                             |
          02_leiden_resolution_sweep.py                                            step03_cluster_and_plot.py
          (binary-searches resolution to hit                                       (fixed resolution -- the value
           target k=6 domain_true / k=8                                             02 just found -- UMAP colored
           cell_type_true; no plots)                                                by ground truth/cluster_label,
                    |                                                                + contingency heatmap)
                    v
          composition_recovery.py
          (scores 02's output vs. soft
           ground-truth composition +
           slice_id leakage; used only
           by the BANKSY/domain-sweep
           scripts, not the plain
           cell-type panels)
                    |
                    v
          ari_recovery_summary/plot_ari_recovery.py
          ari_recovery_summary/plot_domain_vs_celltype.py
          strong_domain_mix/batch_sigma_slide/plot_batch_sigma_slide_final.py
```

## Cell-type clustering and deconvolution

Inputs: `sim_paper/data/figure_2/smaller_sphere/data/<tag>/`.
All clustering panels use the QC-filtered h5ad, including cell.

| Technology | Figure 2 tag |
|---|---|
| cell | `log_mu_-2.5_theta_0.40_jitter0.15_bsigma15` |
| bin8um | `packing_pf0p04_log_mu_0.0_bsigma08` |
| bin16um | `packing_pf0p04_bin16um_log_mu_-2.5_bsigma07` |
| spot | `packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03` |

`weak_domain_mix/pca_harmony/{cell,bin,bin16um,spot}_celltype_panel.sh`
run PCA → Harmony → resolution-matched Leiden against `cell_type_true`.
Results live under `data/figure_3/weak_domain_mix/pca_harmony/<modality>/bs<sigma>/`
(`bin` denotes bin8um here). Existing PCA outputs are reused; changing an
input tag alone does not invalidate an existing result.

`spatial_deconvolution/weak_domain/run_rctd_spot.sh` runs RCTD with the current
Figure 2 QC cell data as the labelled reference and QC spot data as the query
(Figure 2's weak/realistic domain mix); the independent-seed reference scripts
live there too. `spatial_deconvolution/strong_domain/run_rctd_spot_strong_domain.sh`
runs the same RCTD on the strong-domain-mix data (spot bs0.3 query; cell bs1.5
references at seeds 2025/101/202, all from `strong_domain_mix/batch_sigma_slide/`), so it
needs the strong-mix and seed-replication jobs done first.
`average_rctd_seeds.py --config {weak_mix,strong_mix}` averages each
configuration's 3 references. Results: `data/figure_3/spatial_deconvolution/RCTD/`
`weak_mix{,_seed*,_seed_avg}/` and `strong_mix{,_seed101,_seed202,_seed_avg}/`.
Plotting (`plot_rctd_results.py`, `average_rctd_seeds.py`) stays at the top of
`spatial_deconvolution/`.

## Domain recovery and batch-effect sweeps

The current baseline generators are:

```bash
sbatch sim_paper/code/clustering/strong_domain_mix/generate/generate_strong_mix_cell.sh
sbatch sim_paper/code/clustering/strong_domain_mix/generate/generate_strong_mix_bin16um.sh
sbatch sim_paper/code/clustering/strong_domain_mix/generate/generate_strong_mix_spot.sh
```

Each runs generation → QC → BANKSY + Harmony → matched Leiden/ARI →
composition_recovery (including slice-batch leakage diagnostics).
All three use sphere radius 2050 µm. Cell uses 24,207 cells and the current
Figure 2 dispersion; spot also uses its current Figure 2 dispersion.

| Technology | BANKSY lambda | k_geom | Canonical batch_sigma |
|---|---:|---:|---:|
| cell | 0.5 | 60 | 1.5 |
| bin16um | 0.5 | 100 | 0.7 |
| spot | 0.1 | 8 | 0.3 |

Generators write raw/QC inputs under `data/noisy/<technology>_strong_mix_bs<value>/`.
BANKSY matrices, ARI results and scores go under
`data/figure_3/strong_domain_mix/batch_sigma_slide/`, with per-point folders
`<cell|bin16um|spot>/bs<value>/`. Some older baseline raw/QC files were moved
into those result folders; the generators do not automatically reuse those
moved files. Do not assume rerunning a generator uses the same on-disk input
as a downstream consumer without checking its path.

Cell baseline grid: 0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5.
Bin16um generator grid: 0, 0.25, 0.30, 0.35, 0.40, 0.45, 0.7;
its existing 0.05 baseline is also included in seed replication and plotting.
Spot baseline grid: 0, 0.05, 0.1, 0.12, 0.13, 0.14, 0.15, 0.18, 0.2,
0.22, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3; seed work extends it to 0.4 and 0.6.

`../misc/clustering/misc/banksy_batch_compare/` and `../misc/clustering/misc/legacy/` contain historical workflows;
they are not prerequisites for the current baseline generators or plots.
The strong-mix simulations are also intended for STAGATE. Its path audit
is deferred: cell strong-mix paths currently expect raw/QC files in the
Figure 3 results tree, while regenerated cell files are under `data/noisy/`.
No STAGATE changes are part of this update.

## Seed replication and plotting

Baseline seed is 2025; additional simulation seeds are 101 and 202.
These change the simulation draw, not just the clustering initialization.
The seed script matches the baseline generation/BANKSY settings above.

```bash
python sim_paper/code/clustering/strong_domain_mix/batch_sigma_slide/gen_batch_sigma_slide_seed_tasks.py
# Full array for a fresh run (do not duplicate an already submitted array):
sbatch sim_paper/code/clustering/strong_domain_mix/batch_sigma_slide/batch_sigma_slide_seeds.sh
```

The generated TSV contains 72 tasks: cell 18, bin16um 16, spot 38.
Indices 0–67 preserve the original task mapping. Indices 68–71 append cell
0.05 seed101/seed202 and bin16um 0.05 seed101/seed202. Spot 0.4/0.6 include
baseline generation as well as the two extra seeds. Seed simulation files
live under `data/noisy/<cell|bin|spot>_batch_slide_bs<value>_seed<seed>/`;
results live under `data/figure_3/strong_domain_mix/batch_sigma_slide/<cell|bin16um|spot>/bs<value>_seed<seed>/`.

```bash
python sim_paper/code/clustering/strong_domain_mix/batch_sigma_slide/plot_batch_sigma_slide_final.py
```

The plot reads domain ARI JSON files from baseline and seed folders,
groups numerically equivalent batch values (e.g. 0.30 and 0.3), and plots
the available-seed mean ± sample standard deviation. Points with fewer
than three completed seeds are labelled `n=1` or `n=2`; a single seed has
no SD bar. Duplicate entries for the same batch value and seed raise an
error rather than being counted twice. Missing composition-score files do
not exclude a completed ARI result.

Outputs in `data/figure_3/strong_domain_mix/batch_sigma_slide/`:
- `batch_sigma_slide_domain_ari_final.png`
- `batch_sigma_slide_domain_ari_final.csv`: mean, SD, count, included seeds,
  and missing seeds for each technology/batch value.

Regenerate the plot after jobs complete; an interim plot is not a complete
three-seed result. The original array is job `35788619` (68 tasks; pending
at the last audit). Additional control tasks are submitted separately to
avoid duplicating that work; see the run log below.

## Summary plots and environments

`ari_recovery_summary/plot_ari_recovery.py` summarizes Figure 2 cell-type
recovery. `ari_recovery_summary/plot_domain_vs_celltype.py` compares strong-mix
baseline domain ARI at batch_sigma=0 against canonical Figure 2 cell-type
ARI for cell, bin16um and spot. That comparison is separate from the
seed-aggregated batch-effect curve.

Scripts source `_env.sh` (albis-tutorial) or `_env_banksy.sh`
(sim-app-banksy). BANKSY has its own environment because its dependencies
differ from the main analysis stack. RCTD uses `sim_paper/env/rctd`.

## Steps, if you need to run them individually

1. **Generate** the sim dataset: `sim_paper/code/data/generate_simulation_noisy.py`
   (not in this directory -- see `sim_paper/code/data/`).
2. **`00_qc_filter.py`** -- drops fully-empty and near-empty observations
   before clustering (targets neighbor-graph fragmentation from empty
   bins/spots, not sim-vs-real fairness; see the script's docstring).
   ```bash
   python sim_paper/code/clustering/00_qc_filter.py --modality bin --packing-tag <tag>
   ```
3. **Cluster** -- one of:
   - `01.2_pca_harmony.py`: plain PCA -> Harmony (on `obs["slice_id"]`, the batch
     key albis's synthetic batch effects are applied per). Writes PCA and
     before/after-Harmony UMAP diagnostic plots plus PC-pairs plots
     (`pc_pairs.py`, imported, not a standalone CLI).
     ```bash
     python sim_paper/code/clustering/01.2_pca_harmony.py --modality bin --input <qc.h5ad> --output <out.h5ad>
     ```
   - `01_build_banksy_matrix.py` (needs the `sim-app-banksy` env): BANKSY ->
     PCA -> Harmony, staggering each slice's coordinates first (all slices
     share one local x/y range, so an unstaggered neighbor graph would treat
     different slices as spatial neighbors). Writes `obsm["X_pca_harmony"]`
     in the same shape `01.2_pca_harmony.py` does, so
     `step03_cluster_and_plot.py` needs zero changes to consume either.
     ```bash
     python sim_paper/code/clustering/01_build_banksy_matrix.py --modality cell --packing-tag <tag>
     ```
4. **`02_leiden_resolution_sweep.py`** -- Figure 3B: binary-searches the Leiden
   resolution until the cluster count matches the true category count (so
   ARI isn't confounded by over/under-clustering), separately for
   `domain_true` and `cell_type_true`. Runs BEFORE step 5 below, not after --
   see the note under it.
   ```bash
   python sim_paper/code/clustering/02_leiden_resolution_sweep.py --modality bin --packing-tag <tag>
   ```
5. **`step03_cluster_and_plot.py`** -- Leiden or Louvain clustering
   (`--algorithm`, default `leiden`) + UMAP from the pca_harmony pipeline's
   output (`--pipeline` kept as a single-choice flag, `pca_harmony`, since a
   second gene-space Harmony pipeline lived here before being deleted
   2026-08-27; its dead branch was stripped 2026-09-16), plots UMAP
   colored by `domain_true`/`cell_type_true`/`cluster_label`, and (if both a
   `--cluster-key` and a `*_true` ground-truth column are present) a
   ground-truth-vs-predicted UMAP + contingency heatmap.
   ```bash
   python sim_paper/code/clustering/step03_cluster_and_plot.py --modality bin --resolution 0.5
   ```
   Renamed 2026-09-17 from `clustering_plots.py` (no numeric prefix at all)
   to `step03_cluster_and_plot.py` -- `step03` rather than a bare `03_`
   because it's imported as a library (`from step03_cluster_and_plot import
   plot_umap_true_vs_predicted, ...`) by 3 other scripts, and Python can't
   `from 03_x import y` -- digits can't start an identifier, but `step03...`
   is a valid one and still sorts/reads in step order. Every panel script
   actually runs step 4 first to find the right resolution, extracts it from
   `ari_summary_<mod>.json`, then passes that exact value to this step's
   `--resolution` -- a fixed guessed resolution can drift off the true
   category count (found 2026-08-25 on spot: res=0.5 gave 7 clusters vs. 8
   true types) and make this step's qualitative plots look like a mismatch
   that isn't really there.
6. **`plot_ari_recovery.py`** -- grouped bar chart of `cell_type_true` ARI
   across modality, from `weak_domain_mix/pca_harmony/`. Trimmed 2026-09-16 to
   drop the `domain_true`/BANKSY half (it only ever read
   `data/clustering_<tag>/<modality>/banksy_ari_recovery/`, populated
   exclusively by `legacy/*_domain_panel.sh`, last built 2026-08-25) -- for
   domain-vs-celltype comparison see `ari_recovery_summary/plot_domain_vs_celltype.py`
   instead, which reads current data.
   ```bash
   python sim_paper/code/clustering/ari_recovery_summary/plot_ari_recovery.py
   ```

Corresponding `run_*.sh` SLURM wrappers for steps 2-6 (`0-2` array over
spot/bin/cell where relevant) live in `../misc/clustering/misc/` (`run_qc_filter.sh`,
`run_pca_harmony.sh`, `run_clustering_leiden.sh`, `run_ari_vs_ground_truth.sh`)
-- moved from the top level 2026-09-16: nothing in the current pipeline calls
them (every panel/sweep script invokes the `.py` steps directly), and their
only `figure.md` references are old one-off historical reproduce recipes
(Figure 2 write-up, the 2026-09-06 BANKSY tuning sweep), not live
instructions. Still useful if you want to hand-run a custom one-off sweep the
way that work did -- see figure.md's "Reproduce" code blocks for worked
examples.


## Run log

- 2026-09-18: preserved original array `35788619` (indices 0–67). Submitted
  additional array `35789455` (indices 68–71): cell and bin16um at
  batch_sigma=0.05, seeds 101 and 202.
- Verified original TSV rows are unchanged, all 72 tasks are unique, and
  numeric batch grouping, mean/sample SD and incomplete-seed handling work.
  Rendered the updated plot and CSV on existing baseline results; all
  currently available points have n=1. Rerun plotting after the arrays finish.

- 2026-09-18: standardized zero-batch definition to post-resampling X at
  batch_sigma=0 (all systematic multipliers equal one, Poisson variation
  retained). Cell/spot and seed replicates already follow this definition.
  Archived bin16um bs0 and its score under
  `bin16um/bs0_pre_poisson_baseline_20260918/`; regeneration + BANKSY job
  `35789683_0` replaces its historical pre-batch-layer baseline.

## Count-matrix policy (2026-09-18)

Active BANKSY and Figure 4 analyses start from the input h5ad `.X`,
including sigma=0 with Poisson resampling. BANKSY and STAGATE no longer
accept `--use-pre-batch`; old oracle wrappers under `../misc/clustering/misc/` are historical
and will fail on that obsolete flag rather than silently switching layers.
Pre-batch layers remain stored as simulation provenance, not analysis inputs.
The Figure 4 plotting launcher now selects the current strong-mix seed grid
and the main panel excludes pre-batch/oracle metrics and historical runs.
The STAGATE panel now displays means ± sample SD across seeds 101, 202,
and 2025, with individual seed markers and actual batch sigma on the x-axis.
Archived standalone plotting scripts remain historical.

## Current manuscript plotting audit (2026-09-19)

The STAGATE input audit described above as deferred has since been completed;
see `../applications_albis/spatial_clustering/README.md`. The Figure 3 batch
summary CSV now provides the seed coverage for each plotted condition; historical
n=1 run-log entries are snapshots, not current status.

Selected RCTD panels use `spatial_deconvolution/RCTD/weak_mix_seed_avg/` (strong-mix counterpart: `strong_mix_seed_avg/`): estimated
fractions averaged across reference seeds 2025 (canonical), 999999, 314159 before computing
r and RMSE. Spatial zoom displays four types in section 5; the error distribution
covers all eight types. Violin inputs are capped at the pooled 99th percentile
for display; metric calculations and boxplot inputs use the uncapped errors.
The caption should state this display cap if the violin panel is retained.
