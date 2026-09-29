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

**2026-09-23 spot decision:** we dropped the experimental
`--domain-size-factors 0.35 0.6 1.0 1.0 1.6 2.8`. Despite improved
total-count distribution fit, these factors made spatial domains partly
identifiable from library depth. The adopted Figure 2 spot configuration
uses `--base-gene-lognormal -2.25 1.0 --theta 0.25 --theta-jitter 0.10
--batch-sigma 0.3`, with no domain size factors. Figure 3 uses the same
expression settings, including in its strong-mix simulations. Earlier
`dsf*` results are experimental history; see the
[final decision](../../data/figure_2/smaller_sphere/test/README.md).

Eight modality-vs-real pairings, **one script each** under this directory
(`<pairing>.sh` -- the `_with_batch.sh` / `_with_batch_tuned.sh` variants
were retired 2026-08-27 once slice-5 + post-batch became the defaults, and
`spot_vs_breast_cancer_visium.sh` was retired 2026-09-05 when spot's real
reference moved to the two probe panels, see "Superseded output" below):

| Modality | Real reference | Sim config tag |
|---|---|---|
| cell | Xenium `non_diseased_lung` | `log_mu_-2.5_theta_0.40_jitter0.15_bsigma15` |
| cell | Xenium `lung_cancer` | `log_mu_-2.5_theta_0.40_jitter0.15_bsigma15` (same) |
| bin (8um) | Visium HD `breast_cancer_visium_hd` | `packing_pf0p04_bin8um_from16umcfg_log_mu_-2.5_jitter0.6_bsigma07` (2026-09-23: bin16um tissue at 8um) |
| bin (8um) | Visium HD `human_pancreas_visium_hd` | `packing_pf0p04_bin8um_from16umcfg_log_mu_-2.5_jitter0.6_bsigma07` (same) |
| bin (16um) | Visium HD `breast_cancer_visium_hd_16um` | `packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma07` |
| bin (16um) | Visium HD `human_pancreas_visium_hd_16um` | `packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma07` (same) |
| spot | Visium `lymph_node_visium` (probe, CytAssist FFPE) | `packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_bsigma03` |
| spot | Visium `tonsil_visium` (probe, CytAssist FFPE) | `packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_bsigma03` (same) |

**Spot (2026-09-05):** `breast_cancer_visium` (36k whole-transcriptome,
fresh-frozen) was dropped as spot's tuning target on 2026-09-04 and replaced
by the two 18k CytAssist probe references, tuned jointly against both. The
config moved `log_mu` -2.5 -> -2.0 and `theta` 2.0 -> 0.25 (`sphere_r_um`
2050 and `batch_sigma` 0.3 held). Picked by the probe sweeps
(`code/misc/data/misc/sweep_spot_logmu_probe.sh`, `sweep_spot_theta_probe.sh`,
`sweep_spot_logmu_theta_joint_probe.sh`, tabulated by
`summary_spot_logmu_theta_joint.sh`) on the lowest composite = sum of
`|ln(sim/real ratio)|` over {`theta_hat`, `total_counts_median`,
`matrix_zero_frac`}, summed across BOTH probe refs. lymph_node is the strong
half of the fit; tonsil's real HVG-matched `theta_hat` (~0.34) is low enough
that every swept config overshoots its dispersion by >=1.8x.

**Spot (2026-09-06) -- `--theta-jitter` bug fix.** The shared spot config was
generating with `--theta 0.25` and **no `--theta-jitter`**, so
`generate_simulation_noisy.py` fell back to its module default
`NOISY_THETA_JITTER = 1.0`. Per-gene NB dispersion is drawn `~N(theta,
theta_jitter)` = `~N(0.25, 1.0)`, so ~40% of genes drew a value <= 0 and were
floored to `1e-3` (extreme overdispersion). Result: `mean_variance_compare.png`
and `mean_dropout_compare.png` showed the sim genes as **two disjoint clouds** --
a dense upper cloud (the floored genes) far above the NB fit, plus the main
cloud on it. This is the **same artifact fixed for `cell` on 2026-08-25**
(`--theta-jitter 0.15`, cell tag `..._jitter0.15_...`); the fix was never
carried over when spot was retuned onto the probe refs. `--theta-jitter` was
bracketed (0.10 / 0.15 / 0.25, `code/misc/data/misc/sweep_spot_jitter_probe.sh`
job 35536033, tabulated by `summary_spot_jitter.sh`): **0.10** is the largest
value that fully collapses the two clouds into one continuous locus (0.4% of
genes floored vs 40%; per-gene `theta_hat` median ~0.20, smooth unimodal).
0.15 leaves a faint residual cloud, 0.25 clearly reforms it. `theta` itself
is unchanged at 0.25 -- the 2026-09-05 joint `theta` sweep ran with
`jitter=1.0`, and the floored subpopulation drags the median `theta_hat`
down, so `theta` should be re-checked now that jitter is sane (open, see
`figure.md` 2026-09-06). The jitter sweep's `spot_jitter_sweep_0.10` build
(identical generate flags + `--theta-jitter 0.10`) was staged into
`data/figure_2/packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03/`,
no regeneration.

Each script is self-contained: generates the sim data if missing, QC-filters
it (`code/clustering/00_qc_filter.py` -- drops fully-empty and near-empty
observations, a graph-fragmentation fix borrowed from the Figure 3
clustering pipeline, not a real-data-parity correction), then runs
`count_distribution.py` in four modes (see below). Bin/16um real references
are resampled from the same raw 10x `binned_outputs/` at 16um instead of
8um (`code/misc/real_data_qc/misc/run_visium_hd_qc_16um.sh`) -- both resolutions
are canonical Visium HD configurations, not a replacement of one by the
other.

**The same sim `.h5ad` files are used for Figure 3's cell-typing panels**
(`code/clustering/<modality>_celltype_panel.sh`, plus a
`bin16um_celltype_panel.sh`), so Figure 2 and Figure 3 report on
byte-identical data -- the two spot pairings above share one sim file, as do
the two cell / two bin8 / two bin16 pairings, so it is four distinct sim
configs across the eight pairings. Figure 2's scripts and the celltype-panel
scripts each carry a matching generate-if-missing fallback with the identical
config. `spot_celltype_panel.sh` was moved to the new spot tag on 2026-09-05
alongside this retune and resubmitted; its pre-retune outputs are archived at
`data/figure_3/pca_harmony_single_cell/spot_pre_probe_retune_20260905/`, and
the pre-`--theta-jitter`-fix outputs (2026-09-06) at
`.../spot_pre_jitter_fix_20260906/`. The
Figure 4/5 application scripts (`code/applications_albis/**`,
`code/clustering/domain_celltype_composition.py`, the `_strongmix`
derivative) still read the old `packing_pf0p04_log_mu_-2.5_bsigma03` tag,
which stays on disk -- migrating those is a separate pass.

## Three methodological choices, now baked in as defaults

**1. Single representative slice, not pooled.** Every sim dataset has 10
z-plane slices through the simulated tissue sphere (`slice_id` 0-9). Real
reference data is a single tissue section. Pooling all 10 slices is not
apples-to-apples, and empirically two of them (`slice_id` 0 and 9, nearest
the sphere's poles) are geometrically degenerate -- their cross-sectional
disc is much smaller than at the equator, so most observations land
off-tissue (~21% empty bins vs. <0.5% for the 8 interior slices on the
bin16um config). Pooling silently dragged the aggregate stats down.
**Decision: restrict to one central slice**, applied to the sim/primary side
only (real data has no `slice_id`). This is the `count_distribution.py`
**default**; pass `--all-slices` to pool. `slice_id` 0 is the lowest z.
**2026-09-29: changed from `slice_id=5` (5th from the top) to `slice_id=4`
(5th from the bottom)**, so Figure 2 and Figure 5 (5A/5B/5D, Supplementary
Figure 6) use the same physical slice. With 10 slices there is no single
middle slice; 4 and 5 are its two mirror-image central sections.
Slice-5 outputs are archived at
`data/figure_2/misc/archive/smaller_sphere_plots_slice_id5_20260929/`.

**2. Batch effect: included (post-batch counts).** The sim's synthetic
per-slice technical batch shift (`batch_sigma`) is included in the
count-distribution stats -- `count_distribution.py` uses `X` (post-batch)
uniformly by **default**; pass `--no-batch-effect` to fall back to the
`counts_pre_batch` layer. Real data has no "pre-batch" version -- whatever
technical noise a real section carries is just baked into its counts -- so
including the sim's batch effect is the fairer comparison now that we're
down to one slice. It is *not* a small effect (bin16um `theta_hat` at slice
5: ~1.2 pre-batch vs. ~0.18 post-batch).

**`batch_sigma` is set per modality** (cell 1.5, bin 8um 0.8 -> 0.7 on 2026-09-23 when it adopted the bin16um tissue, bin 16um 0.7,
spot 0.3), finalized 2026-08-27. The value is chosen on Figure 3's
pre/post-Harmony slice-separation demo -- each modality's aggregation
footprint dilutes the fixed-size per-slice shift differently, so one global
value would be invisible for some modalities and overwhelming for others --
then confirmed here to still match the real count distribution (spot's and
cell's `theta_hat` in particular are sensitive to it). Sweep scripts:
`code/misc/clustering/misc/batch_sigma_sweep.sh` and `batch_sigma_opt_round2.sh`.

**3. Real platform capture window ("realwindow"), baked in 2026-08-31.**
`generate_simulation_noisy.py` no longer scales the capture window with
`--sphere-r-um`. With `--capture-window-um` unset it uses the real instrument
window: 6.5 x 6.5 mm for bin/spot (Visium / Visium HD), 12 x 24 mm for cell
(Xenium, via `xenium_capture_window_um`). The ~2050 um simulated tissue disc
now sits inside the full Visium(HD) window with a wide empty border, so
**~77% of bins/spots are off-tissue before QC**. Off-tissue observations get
`domain_true` / `cell_type_true` = `"unassigned"` and `obs["is_empty"] = True`
(previously `np.argmax` of an all-zero fraction row silently labelled them
`D0`, contaminating that domain's ground truth and producing the polar-slice
border artifact). `00_qc_filter.py` drops them, so the QC'd sim is a clean
6-domain / 8-cell-type dataset. This mirrors real Visium HD, which ships the
whole capture window and excludes off-tissue bins via `in_tissue` + QC. cell
is unaffected in practice (cells are only placed within the tissue sphere,
r = 6000 < 12 mm, so there are no off-tissue cells). Pre-realwindow sim data
is archived at `data/figure_2_oldwindow_20260831/`; the prior comparison
outputs at `data/count_distribution/sweeps/figure_2_oldwindow_20260831/`.

## The four comparison modes

Each pairing script runs `count_distribution.py` four times, into
`<OUT_ROOT>/{full_panel,hvg_matched,qc_filtered,qc_and_hvg_matched}/`
under `sim_paper/data/count_distribution/figure_2/<modality>_vs_<real>/`:

| Mode | Sim input | `--match-panel-size` | Use for |
|---|---|---|---|
| `full_panel` | raw (not QC'd) | no | `mean_variance`/`mean_dropout` dispersion fit -- full gene panel, unfiltered; HVG-selection or QC would bias the dispersion estimate itself. |
| `hvg_matched` | raw (not QC'd) | yes | `sparsity_summary`/empty-row-rate -- panel-matched but not QC'd, since QC-filtering sim here would hide the zero-inflation these stats exist to report. |
| `qc_filtered` | QC'd | no | `total_counts` (secondary use) with QC'd sim, full gene panel. |
| `qc_and_hvg_matched` | QC'd | yes | `total_counts`/`genes_per_cell` -- both corrections together, the fair like-for-like for these two panels specifically. |

**Since the realwindow switch (2026-08-31), `full_panel` and `hvg_matched`
for bin and spot are dominated by off-tissue empty rows** (sim
`total_counts_median` = 0, `empty_row_frac` ~0.69 at slice 5) and are not
used in the figure -- the bin/spot panels are read only from `qc_filtered` /
`qc_and_hvg_matched`. The two raw modes stay meaningful for **cell** (no
off-tissue observations) and are still generated for all modalities as
diagnostics.

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
sbatch sim_paper/code/count_distribution/<pairing>.sh
```

Each is idempotent for the generate/QC steps (skips if the sim `.h5ad`
already exists) but always reruns all 4 `count_distribution.py` modes.
To change the sim config, edit `SIM_TAG` + the `generate_simulation_noisy.py`
flags together at the top of the relevant script(s) and resubmit -- doesn't
touch the other pairings. Keep the matching
`code/clustering/<modality>_celltype_panel.sh` in sync.

`code/misc/data/misc/generate_figure2_batch_tuned.sh` is a convenience SLURM
array that (re)generates all four sim configs at once into
`data/figure_2/<tag>/`; the per-pairing scripts don't need it but it's the
fastest way to rebuild everything.

For a quick one-off check at a different slice or without the batch effect,
pass `--all-slices` or `--no-batch-effect` to `count_distribution.py`.

## Superseded output (archived)

Retired 2026-08-27 into `sim_paper/data/count_distribution/sweeps/`:
- `figure_2_prebatch_slice5/` -- slice 5, pre-batch (the old default).
- `figure_2_with_batch_flat022/` -- slice 5, post-batch, flat `batch_sigma=0.22`.
- `figure_2_with_batch_tuned_presweep/` -- slice 5, post-batch, the earlier
  0.22/0.5/0.5/0.15 per-modality values before the 2026-08-27 re-sweep.

Retired 2026-08-31 (realwindow switch):
- `sim_paper/data/figure_2_oldwindow_20260831/` -- the 4 sim `.h5ad` tags as
  generated under the old sphere-scaled tight capture window.
- `sim_paper/data/count_distribution/sweeps/figure_2_oldwindow_20260831/` --
  all 7 comparison outputs from those tags.
- `code/clustering/*_celltype_panel.sh` (Figure 3) still point at the
  now-archived old-window tags and need the same regeneration + rerun --
  tracked in `figure.md` / `DATA_VERSIONS.md`, not done in this pass.

Retired 2026-09-05 (spot moved to probe references):
- `code/count_distribution/spot_vs_breast_cancer_visium.sh` -- **deleted.** Its
  outputs stay on disk at
  `data/count_distribution/figure_2/spot_vs_breast_cancer_visium/` (last run
  2026-09-01, the old `packing_pf0p04_log_mu_-2.5_bsigma03` config).
- `data/figure_2/packing_pf0p04_log_mu_-2.5_bsigma03/` -- **kept on disk**, not
  archived: still read by the Figure 4/5 application scripts. The new spot tag
  `packing_pf0p04_log_mu_-2.0_theta_0.25_bsigma03` sits alongside it, staged
  from the `spot_joint_logmu_-2.0_theta_0.25` sweep build.
- `data/figure_3/pca_harmony_single_cell/spot_pre_probe_retune_20260905/` --
  the pre-retune Figure 3 spot cell-typing outputs, moved aside so
  `spot_celltype_panel.sh` (not tag-scoped at `CLUSTER_ROOT`) recomputes on
  resubmit instead of skipping.

Retired 2026-09-06 (spot `--theta-jitter` bug fix, see "Spot (2026-09-06)"
above):
- `data/count_distribution/sweeps/figure_2_spot_pre_jitter_fix_20260906/`
  `spot_vs_{lymph_node,tonsil}_visium/` -- the pre-fix Figure 2 spot
  comparison outputs (two disjoint gene clouds in the mean_variance /
  mean_dropout panels).
- `data/figure_3/pca_harmony_single_cell/spot_pre_jitter_fix_20260906/` --
  pre-fix Figure 3 spot cell-typing.
- `data/figure_4/alignment/STAIR/spot_pre_jitter_fix_20260906/` -- pre-fix
  Figure 4C STAIR spot alignment.
- `data/figure_2/packing_pf0p04_log_mu_-2.0_theta_0.25_bsigma03/` -- the
  pre-fix (`theta_jitter=1.0`) spot sim data, kept on disk (still read by the
  un-migrated Figure 4/5 application scripts, same as the older
  `..._log_mu_-2.5_bsigma03` tag).

## Known open gaps (not addressed by this methodology, tracked in figure.md)

- Pancreas is a known-bad calibration target for bin (extreme real
  overdispersion from a few dominant hormone genes) -- not expected to
  match well regardless of slice/batch choice.
- Spot's `genes_per_cell` still mismatches real even panel-matched --
  open issue, unrelated to slice/batch. After the 2026-09-05 probe retune,
  qc_and_hvg_matched `genes_per_cell` median is ~274 sim vs ~394 lymph_node
  (undershoots) but ~258 tonsil (close) -- not in the composite the sweep
  optimizes.
- **RESOLVED 2026-09-13:** Spot's `theta=0.25` re-checked with `theta_jitter=0.10`
  fixed (the re-bracket this item asked for). A lower `theta=0.15` scores
  better on the composite (1.773 vs 2.764) but, holding `jitter=0.10` fixed,
  reintroduces a small floored-gene upper-cloud artifact (0.36% -> 6.8% of
  genes) since `theta`/`theta_jitter` aren't independent -- the jitter would
  need rescaling (~0.06) to match, which wasn't verified. Decision: **kept
  `theta=0.25`/`jitter=0.10`**, prioritizing the already-verified
  artifact-free mean-variance/mean-dropout plots over the marginal composite
  gain. Tonsil `theta_hat` (~0.34 real HVG-matched) remains the limiting term
  regardless of `theta` -- a known ceiling, not something this bracket could
  fix. See `figure.md`, 2026-09-13, and
  `data/count_distribution/sweeps/spot_theta_rebracket_probe/` for the sweep.
  Figure 2 spot config is final.
- QC-parity gap: real bin/spot data goes through SpotSweeper local-outlier
  QC before comparison; sim only gets the minimal `00_qc_filter.py` pass
  (empty/near-empty only). Parked, not addressed here.
- cell was retuned 2026-08-27 (log_mu -2.3, theta 0.40) so the WITH-batch
  theta_hat lands ~0.083 (real Xenium ~0.16) instead of ~0.056, and
  cell_type_true ARI recovers to 0.64. Remaining soft spot: median
  total_counts overshoots (~327 vs real ~90) -- the batch_sigma=1.5 needed
  for the visible Figure 3 Harmony demo inflates it. Accepted for the shared
  dataset; revisit if Figure 2 cell total-counts fidelity becomes priority.
