# Data versions

Tracks the different `sim_paper/data/` simulation runs, what parameters
distinguish them from the manuscript baseline, and which ones are decided
vs. still open. All variants share the base tissue geometry unless noted
(`sphere_R_um=6000`, `n_cells=600_000`, `capture_window_um=(6500,6500)`,
`bin_size_um=8.0`, `n_slices=10`, `n_domains=6`, `n_cell_types=8`).

## Status by modality

| modality | status | decision |
|---|---|---|
| cell | **decided** | `batch_effect_moderate` (`batch_sigma=0.8`) — UMAP looks good, approved 2026-08-20 |
| bin (Visium HD) | **blocked** | UMAPs unusable for publication — see "Known issue" below; fix in trial |
| spot (Visium) | **blocked** | same root cause as bin, less severe — see "Known issue" below. Figure 2 count-distribution config retuned 2026-09-05, see section below |

## 2026-09-06: spot `--theta-jitter` bug fix (two-cloud artifact)

Figure 2's spot `mean_variance_compare.png` / `mean_dropout_compare.png`
(vs both probe refs) showed the sim genes as **two disjoint clouds** -- a
dense upper cloud far above the NB fit plus the main cloud on it. Root
cause: the shared spot config generated with `--theta 0.25` and **no
`--theta-jitter`**, so `generate_simulation_noisy.py` used its module
default `NOISY_THETA_JITTER = 1.0`. Per-gene NB dispersion `~N(0.25, 1.0)`
-> ~40% of genes floored to `1e-3` -> the upper cloud. **Same artifact
fixed for `cell` on 2026-08-25** (`--theta-jitter 0.15`); never propagated
to spot through the 2026-09-05 probe retune.

- **Fix:** `--theta-jitter 0.10`. Bracketed 0.10/0.15/0.25
  (`code/misc/data/misc/sweep_spot_jitter_probe.sh`, job 35536033,
  `summary_spot_jitter.sh`). 0.10 is the largest value that fully collapses
  the two clouds (0.4% genes floored vs 40%; per-gene `theta_hat` median
  ~0.20, smooth). 0.15 leaves a faint residual, 0.25 clearly reforms the
  cloud. The `hvg_matched` composite nominally favours 0.25 but that is the
  floored subpopulation gaming the metric (drags median `theta_hat` toward
  tonsil's low ~0.34), so 0.10 was picked on the "collapse the clouds" goal.
- **New tag:** `packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03`,
  staged from `data/noisy/spot_jitter_sweep_0.10/` (identical generate flags
  + `--theta-jitter 0.10`, rng seeded) into `data/figure_2/<tag>/` -- no
  regeneration. The pre-fix `..._theta_0.25_bsigma03` tag (`theta_jitter=1.0`)
  stays on disk (un-migrated Fig 4/5 scripts).
- **Scripts:** `SIM_TAG` bumped + `--theta-jitter 0.10` added to the generate
  fallback in `code/count_distribution/spot_vs_{lymph_node,tonsil}_visium.sh`
  and `code/clustering/spot_celltype_panel.sh`; `code/applications_albis/
  alignment/3D_stair.py` spot h5ad path bumped. `theta` itself unchanged
  (0.25) -- re-check pending (2026-09-05 joint `theta` sweep ran at
  `jitter=1.0`).
- **Archived:** `data/figure_3/pca_harmony_single_cell/spot_pre_jitter_fix_20260906/`,
  `data/figure_4/alignment/STAIR/spot_pre_jitter_fix_20260906/`,
  `data/count_distribution/sweeps/figure_2_spot_pre_jitter_fix_20260906/{spot_vs_lymph_node_visium,spot_vs_tonsil_visium}/`.
- **Rerun jobs (submitted 2026-09-06):** **35548398** `spot_vs_lymph_node_visium`,
  **35548399** `spot_vs_tonsil_visium` (Fig 2), **35548400** `spot_celltype_panel`
  (Fig 3), **35548401** `3D_stair --dataset spot` (Fig 4C). bin/cell unaffected
  (different modality; their Fig 2/3 panels and the bin16um/cell STAIR runs
  are byte-identical, not rerun). Check `sacct -j 35548398-35548401` +
  `qc_and_hvg_matched/comparison_summary.json` /
  `ari_recovery_qc/ari_summary_spot.json` /
  `figure_4/alignment/STAIR/spot/adata_results/metrics.json`.
- **Fig 4C assembled panel re-run 2026-09-06:** `plot_figure4c.py`
  (`--color-by slice_id / domain_true / cell_type_true`) → refreshed
  `data/figure_4/alignment/plots/figure4c_*` off the new STAIR spot h5ad.
  bin16um / cell columns unchanged.
- **Figure 4A (SVG: nnSVG / scBSP / SPARK-X) and Figure 4B (STAGATE) are
  DELIBERATELY left on the OLD spot tag** `packing_pf0p04_log_mu_-2.5_bsigma03`
  (pre the 2026-09-05 probe retune, pre the jitter fix). **User decision
  2026-09-06: keep as-is, do not migrate.** Rationale: the 4A/4B results are
  about spatial-signal strength / domain-mix coupling, not spot
  count-distribution fidelity, so the spot count tuning doesn't change the
  conclusions; and the `_strongmix` / `_lowbatch` families those panels
  depend on have no `jitter0.10` siblings. **Consequence: Fig 4A/4B spot is
  NOT the same spot dataset as Fig 2 / Fig 3 / Fig 4C** — this is now noted
  in the `DATASETS` block of `3D_nnsvg.py`, `3D_scbsp.py`, `3D_stagate.py`.
  bin16um and cell are on the current shared tags across all of Figure 4.

## 2026-09-05: spot retuned onto the two CytAssist probe references

`breast_cancer_visium` (36k whole-transcriptome, fresh-frozen) was dropped as
spot's Figure 2 tuning target on 2026-09-04 and replaced by two 18k
probe-based CytAssist FFPE references, `lymph_node_visium` +
`tonsil_visium` (real-data QC: `code/real_data_qc/run_visium_{lymph,tonsil}_qc.sh`).
Spot's shared Fig 2 / Fig 3 sim config:

`packing_pf0p04_log_mu_-2.5_bsigma03` → **`packing_pf0p04_log_mu_-2.0_theta_0.25_bsigma03`**
(`--base-gene-lognormal -2.0 0.7 --theta 0.25 --batch-sigma 0.3 --sphere-r-um 2050`).
`log_mu` -2.5→-2.0 and `theta` 2.0→0.25; `sphere_r_um` and `batch_sigma`
held (`batch_sigma` stays fixed — it is set on the Fig 3 Harmony demo).

- **Sweeps** (`code/misc/data/misc/sweep_spot_logmu_probe.sh`,
  `sweep_spot_theta_probe.sh`, `sweep_spot_logmu_theta_joint_probe.sh`;
  tabulated by `summary_spot_logmu_theta_joint.sh`): winner = lowest
  composite = Σ `|ln(sim/real ratio)|` over {`theta_hat`,
  `total_counts_median`, `matrix_zero_frac`}, summed across BOTH probe refs.
  `-2.0/0.25` scored 1.933 vs next-best `-2.0/0.35` at 2.251. hvg_matched:
  lymph_node `theta_hat` 1.14x / `total_counts` 0.90x / zero-frac 1.52x
  (strong); tonsil `theta_hat` 2.52x / `total_counts` 1.30x / zero-frac 1.10x
  (tonsil's real HVG-matched `theta_hat` ~0.34 is low enough that every swept
  config overshoots its dispersion ≥1.8x).
- **Sim data staged**, not regenerated: `data/noisy/spot_joint_logmu_-2.0_theta_0.25/`
  (the joint-sweep build, byte-identical generate flags, rng seeded) moved to
  `data/figure_2/packing_pf0p04_log_mu_-2.0_theta_0.25_bsigma03/`.
- **Scripts:** `code/count_distribution/spot_vs_lymph_node_visium.sh` rewritten
  + new `spot_vs_tonsil_visium.sh`; `spot_vs_breast_cancer_visium.sh` deleted.
  `code/clustering/spot_celltype_panel.sh` moved to the new tag and
  resubmitted (Figure 3 spot cell-typing recomputed). All three jobs
  submitted 2026-09-05.
- **Archived:** `data/figure_3/pca_harmony_single_cell/spot_pre_probe_retune_20260905/`
  (pre-retune Fig 3 spot outputs, moved aside so the not-tag-scoped
  `CLUSTER_ROOT` recomputes on resubmit).
- **Figure 4C STAIR** (`code/applications_albis/alignment/3D_stair.py`) spot
  entry bumped to the new tag (`…/simulation_spot_z_qc.h5ad`) + resubmitted
  all 3 datasets (jobs **35527576** bin16um, **35527577** spot, **35527578**
  cell); pre-retune outputs archived
  `data/figure_4/alignment/STAIR/{spot,bin16um,cell}_pre_probe_retune_20260905/`.
- **NOT migrated:** the other Figure 4/5 application scripts
  (`code/applications_albis/svg/**`, `spatial_clustering/3D_stagate.py`,
  `code/clustering/domain_celltype_composition.py`, the
  `packing_pf0p04_log_mu_-2.5_bsigma03_strongmix` derivative) still read the
  old spot tag, which stays on disk. **2026-09-06 update: for Fig 4A (SVG)
  and Fig 4B (STAGATE) this is now a deliberate keep, not a pending TODO —
  see the 2026-09-06 section above.**

## 2026-08-31: realwindow is now the Figure 2 default

`generate_simulation_noisy.py` no longer scales the capture window with
`--sphere-r-um`. Unset `--capture-window-um` now uses the real instrument
window (`6500x6500` um bin/spot, `12000x24000` um cell via
`xenium_capture_window_um`). For the Figure 2 bin/spot tags
(`sphere_r_um=2050`) the old rule shrank the window to ~2221 um cropped tight
to the tissue disc; the new rule leaves the ~2050 um disc inside the full
6.5 mm window with a wide empty border — **~77% of bins/spots are off-tissue
before QC**. Off-tissue observations get `domain_true`/`cell_type_true` =
`"unassigned"` + `obs["is_empty"]=True` (albis `_labels_from_fracs`
`empty_label`; previously argmax silently made them `D0`) and are dropped by
`00_qc_filter.py`. cell is unaffected (no off-tissue cells).

- **All 7 `code/count_distribution/*_vs_*.sh` adapted** to the self-contained
  generate-into-`data/noisy/` → QC → move-into-`data/figure_2/` pattern, real
  window (no `--capture-window-um`), same SIM_TAG names, all 4 panel modes
  kept. Only `qc_filtered`/`qc_and_hvg_matched` are used in the figure for
  bin/spot; `full_panel`/`hvg_matched` there are now empty-swamped (still
  valid for cell).
- **4 canonical tags regenerated** under the real window:
  `packing_pf0p04_log_mu_0.0_bsigma08` (bin8),
  `packing_pf0p04_bin16um_log_mu_-2.5_bsigma07` (bin16),
  `packing_pf0p04_log_mu_-2.5_bsigma03` (spot),
  `log_mu_-2.3_theta_0.40_jitter0.15_bsigma15` (cell). bin16 + spot reused the
  clean 2026-08-30 realwindow builds (promoted into `data/figure_2/`); bin8 +
  cell regenerated via `code/misc/data/misc/generate_figure2_batch_tuned.sh`
  (array 0-3, idempotent skip for the two already present), the 7 comparisons
  chained `--dependency=afterok` on the matching array task.
- **Archived:** pre-realwindow sim data → `data/figure_2_oldwindow_20260831/`;
  the 7 old comparison outputs →
  `data/count_distribution/sweeps/figure_2_oldwindow_20260831/`. Deleted:
  `bin16um_realwindow_vs_*.sh`, `generate_realwindow_bin_spot.sh`, the
  `_realwindow*` / `_d0bug_prealbisfix` data dirs.
- **Follow-up (NOT done this pass):** `code/clustering/*_celltype_panel.sh`
  (Figure 3) point at the same 4 tags and need the same regeneration + rerun
  so Fig 2 / Fig 3 stay byte-identical.
- **realwindow bin16 `qc_and_hvg_matched` vs real breast_cancer_visium_hd_16um:**
  total_counts median 175 vs 179, genes/cell 71 vs 80, zero-frac 0.867 vs
  0.852, theta_hat 0.17 vs 0.37 (still ~2x under — pre-existing bin16
  dispersion gap, unchanged by realwindow).

## Known issue: bin/spot zero-inflation from sparse tissue packing

All bin/spot-level runs to date (any `theta`, any `batch_sigma`) inherit the
same base geometry, which packs only 600k cells into a `sphere_R_um=6000`
sphere — a 3D packing fraction of ~0.16% vs. real tissue's near-total
packing. A fixed 8µm bin grid laid over that mostly lands in empty
interstitial space:

- bin: 66.3% of bins are fully empty (0 genes detected) vs. 0.0% in real
  `breast_cancer_visium_hd`
- spot: 9.8% of spots are fully empty vs. ~0% in real data

This shows up directly in the UMAPs as a large ring/fringe of structurally
arbitrary, all-zero points completely separate from the real
cell-type-driven clusters (visible in `clustering_moderate_batch/bin` and
`.../spot` `pca_harmony*/plots/umap_pca_harmony_before_after_*.png`). Raising
`theta` (NB dispersion) does not fix this — it was tested up to 200 with
almost no effect (bin `theta_hat` 0.05 → 0.07) — because the problem is
geometric (empty grid cells), not per-molecule count noise.

**In trial (2026-08-20):** shrinking `sphere_R_um`/`capture_window_um` (and
proportionally the other length-scale params: `core_fuzz_width_um`,
`max_shift`) while holding `n_cells=600_000` fixed, to raise 3D packing
fraction. Local 50k-cell sweeps found no single packing fraction hits both
"zero empty bins" and "correct counts/bin" simultaneously — raising density
kills empty bins (0% by ~10% packing) but overshoots real counts/bin by
5-9x, because closer cell spacing makes neighboring cells' molecule clouds
overlap into the same bin (confirmed this is about spacing, not the
placement algorithm: `allow_cell_overlap=True` gave near-identical numbers
at every packing fraction tested). Settled on ~1.5% packing as a "close
enough, not identical" compromise: cuts bin empty-bins from 66.3% to ~7%,
counts/bin within ~1.6x of real (641 vs 405). Local 50k-cell UMAP check at
this setting (no more empty-bin ring/fringe artifact) confirmed the fix
works qualitatively.

`generate_simulation_noisy.py` now takes `--sphere-r-um`,
`--capture-window-um`, `--core-fuzz-width-um`, `--max-shift`, `--n-cells` as
overridable flags (previously hardcoded) so this is a taggable sweep like
`--theta`. Production run (600k cells, full 10-slice geometry,
`sphere_r_um=2863.1`) submitted 2026-08-20 as `--out-tag packing_pf0p015` →
`data/noisy/packing_pf0p015/` (bin + spot only; cell isn't affected by this
issue), SLURM job 35037599, with `pca_harmony` clustering (job 35037630)
chained via `--dependency=afterok`. Output: `data/clustering_packing_fix/`.

**Result (2026-08-20):** `pf0p015` completed. The full-scale run undershot
the ~7% empty-bin prediction from the local 50k-cell sweep -- actual bin
empty fraction is **16.87%** (vs 0.0% real), nonempty mean 521 counts/bin
(1.3x real ~405). The UMAP ring/fringe artifact is still visible at full
scale (`clustering_packing_fix/bin/pca_harmony/plots/umap_pca_harmony_before_after_by_cell_type_true.png`)
-- the local sweep's calibration does not transfer to full-scale geometry,
so further tuning is being done directly at full scale instead of via local
sweeps. Spot fared better: empty fraction 9.8% -> 3.44%.

**Spot border-artifact finding:** the small fringe of stray points in the
spot post-Harmony UMAP (`clustering_packing_fix/spot/pca_harmony/plots/umap_pca_harmony_before_after_by_cell_type_true.png`)
is the *same* empty-space mechanism, concentrated at the poles: 100% of
empty spots are in `slice_id` 0 and 9 (the two polar Z-slices, 17.2% empty
there vs 0% in slices 1-8), and 100% are `domain_true == D0` (D0 makes up
33% of spots in the polar slices vs as little as 5.9% mid-stack). Not a bug
-- polar slices cut a thinner tissue cross-section, so a fixed spot grid
there lands in more empty space. Expected to shrink as packing fraction
rises; no dedicated fix (e.g. trimming polar slices) attempted yet since
`albis` has no exposed z-margin/pole-trim parameter, only `n_slices`.

**In-progress sweep (submitted 2026-08-20 evening, unattended overnight):**
recalibrating packing fraction directly at full scale (bypassing the local
sweep) since it mispredicted both the empty-bin fraction and the
counts/bin overshoot. Three packing fractions, full spot+bin+cell +
`pca_harmony` clustering each, chained via `--dependency=afterok`:

| tag | `sphere_r_um` | target packing | generation jobs | clustering jobs | output |
|---|---|---|---|---|---|
| `packing_pf0p03` | 2258.5 | ~3% | 35046369 (array 0-2) | 35046374/75/76 (spot/bin/cell) | `data/clustering_packing_pf0p03/` |
| `packing_pf0p04` | 2050 | ~4% | 35044702 (bin, done) + 35046370 (spot+cell) | 35046377/78/79 | `data/clustering_packing_pf0p04/` |
| `packing_pf0p06` | 1800 | ~6% | 35044703 (bin) + 35046371 (spot+cell) | 35046380/81/82 | `data/clustering_packing_pf0p06/` |

Bin-only quick check at `pf0p04` (2026-08-20, before the full sweep above
was launched): empty fraction **5.84%**, nonempty mean 898 counts/bin
(2.2x real) -- better than `pf0p015` on emptiness, worse on overshoot, but
still well short of the 5-9x overshoot the old local sweep predicted at
high packing. `pf0p06` and `pf0p03` full results pending.

Next: once all three land, compare empty-bin fraction, counts/bin
overshoot, and the UMAP ring/fringe artifact (bin) and polar-slice fringe
(spot) across all three packing fractions to pick a final value.

## Runs

| dir | generator | theta | theta_jitter | noise_scale | base_gene_lognormal | batch_sigma | cell `theta_hat` vs real (~0.15-0.16) | notes |
|---|---|---|---|---|---|---|---|---|
| `data/simulation_*.h5ad` (top-level) | `code/data/generate_simulation.py` | 25 (default) | 2.0 (default) | 0.9 (default) | (0.7, 0.7) | 0.22 | 1.56 — far from real | manuscript baseline |
| `data/batch_effect_moderate/` | `code/data/generate_simulation_batch_effect_moderate.py` | 25 (default) | 2.0 (default) | 0.9 (default) | (0.7, 0.7) | **0.8** | 1.56 | **cell: approved.** bin/spot: blocked, see known issue |
| `data/batch_effect_high/` | `code/data/generate_simulation_batch_effect_high.py` (renamed 2026-08-20, was `generate_simulation_batch_effect.py`) | 25 (default) | 2.0 (default) | 0.9 (default) | (0.7, 0.7) | **2.0** | 1.56 | too strong — pushes slices into fully separate clusters, not for publication |
| `data/noisy/` (no tag) | `code/data/generate_simulation_noisy.py` | 2.0 | 1.0 | 1.3 | (0.7, 0.7) | 0.22 | 0.77 | dispersion sweep, first pass |
| `data/noisy/log_mu_-2.5/` | same, `--base-gene-lognormal -2.5 0.7` | 2.0 | 1.0 | 1.3 | (-2.5, 0.7) | 0.22 | 0.76 | lowering baseline expression alone barely moved theta_hat |
| `data/noisy/log_mu_-2.5_theta_0.5/` | same, `--theta 0.5` | 0.5 | 1.0 | 1.3 | (-2.5, 0.7) | 0.22 | 0.27 | |
| `data/noisy/log_mu_-2.5_theta_0.25/` | same, `--theta 0.25` | 0.25 | 1.0 | 1.3 | (-2.5, 0.7) | 0.22 | **0.14** | closest cell-level dispersion match to real Xenium lung to date |
| `data/noisy/theta_200/` (bin only) | same, `--theta 200` | 200 | 1.0 | 1.3 | (0.7, 0.7) | 0.22 | n/a (bin: 0.07) | tested whether raising theta alone fixes bin overdispersion — it doesn't (see known issue) |

## Open gaps

- The `noisy` (dispersion) and `batch_effect` (batch strength) sweeps are
  still separate, un-merged experiment tracks: no single dataset yet combines
  the best dispersion fit (`log_mu_-2.5_theta_0.25`, cell `theta_hat`=0.14)
  with the approved batch strength (`batch_sigma=0.8`).
- No decision yet on which real dataset to benchmark bin/spot against post-fix
  (real refs on hand: `breast_cancer_visium_hd`, `human_pancreas_visium_hd`
  for bin; `non_diseased_lung`, `lung_cancer` for cell).

## Planned: consolidate the generator scripts (target week of 2026-08-24)

Provenance is currently split across 5 scripts with different, overlapping
CLI surfaces -- `generate_simulation.py` (manuscript baseline, no flags),
`generate_simulation_batch_effect_high.py` / `_moderate.py` (hardcoded
radius/geometry, no override), `generate_simulation_noisy.py` (full
dispersion + geometry flags), `generate_simulation_packing.py` (geometry
flags only, dispersion hardcoded). This is why `batch_effect_high` had no
`--sphere-r-um` and why `packing.py` was forked from `noisy.py` tonight
(2026-08-20) instead of just adding a flag to an existing script.

For publication readiness, plan to:
1. Merge into one script exposing every parameter (theta/theta_jitter/
   noise_scale/batch_sigma/base_gene_lognormal/sphere_r_um/n_cells/
   cell_radius_kwargs/...) as a flag, each defaulting to the manuscript
   baseline value. The current named variants become documented flag
   presets/examples rather than separate files with different hardcoded
   subsets.
2. Have the script print its fully-resolved parameter set (already built
   for `uns['sim_params']`) at the top of every run, not just the flags
   that were explicitly passed -- so every SLURM log is a complete,
   self-contained record of what ran, regardless of what was left at
   default.
