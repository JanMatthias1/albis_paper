# Dedicated independent-offset cross-modality experiment

## Native ALBIS reproduction of the stronger-offset condition

All new workflow code lives in this directory; the ALBIS package is used without
modification. `generate_native_offsets.py` calls
`albis.simulate_3d_molecule_sphere_multires` with explicit manuscript presets and
native `max_shift=3075`, `max_deg=270`, `base_seed_unaligned=12345`,
`sync_unaligned_seed=False`, and tissue/expression `seed=2025`.

The shift bound applies independently to each XY translation component. With
the same perturbation seed, changing 1025 to 3075 gives the same rotation and
three times the translation. It does **not** draw a new direction. Change
`--base-seed-unaligned` for new offset directions without changing the tissue;
change `--seed` to change the tissue/expression realization. Neither a single
post-hoc condition nor a single native condition establishes general robustness.

`validate_native_offset_equivalence.py` replayed native draws on the existing
truth for all three modalities of section 5. Maximum coordinate differences
from the saved shift3x inputs were 0.000365 µm or smaller (float32/fit rounding).
The report is `native_offsets_3075_seed12345/native_equivalence.json`. This is
a coordinate-equivalence check, not evidence that a full regeneration has run.

From the project root, submit full generation of all ten sections for each
modality, followed by section-5 STAIR, the approved plots, reference RMSE, and
the two-pair averaged gene metrics:

```bash
bash sim_paper/code/applications_albis/cross_modality_alignment/run_native_offsets.sh submit
```

The default destination is `data/figure_4/cross_modality_alignment/native_offsets_3075_seed12345/`.
The existing `independent_offsets_shift3x/` remains the historical result.
The runner accepts positional overrides after `submit`: output directory,
max shift, perturbation seed, tissue seed. Generation refuses to overwrite a
modality directory, and submission refuses to duplicate a recorded submission.

To inspect settings without writing data:

```bash
sim_paper/env/albis-tutorial/bin/python \
  sim_paper/code/applications_albis/cross_modality_alignment/generate_native_offsets.py \
  --modality bin16um --outdir /tmp/unused --print-config
```

The native recipe writes raw/QC h5ad files, full simulator arguments, source
hashes, software versions and generation status per modality. The runner saves
SLURM logs/job IDs and a source snapshot. The gene evaluator additionally reads
the existing benchmark's `feature_similarity.py` (its path/hash is recorded in
metric provenance). `test_native_offsets.py` checks the parameter/seed behavior;
the existing ALBIS API tests and a small native write/read smoke test also passed.
These additions are local working-tree files until committed to the paper repo;
the recorded ALBIS revision alone does not contain this paper workflow.

### Completed replacement — 2026-09-20

User authorized native regeneration to replace the previous `independent_offsets_shift3x/`
data and plots. Generation array **35816469**, STAIR/plots/metrics
**35816470**, and replacement **35816471** all completed with exit code 0.
All ten sections were regenerated; alignment and evaluation cover section 5.
The previous directory is archived as `independent_offsets_shift3x_before_native_20260920T162259708032Z/`.
The first submission (35816460) failed before generation because
SLURM's spool path was mistaken for the code path; this is fixed. Its blocked
dependent jobs 35816461/35816462 were cancelled and its logs/source snapshot retained.
The replacement only runs after
both upstream steps succeed and verifies generation manifests and required outputs.
`promote_native_offsets.py` archives the old directory under a dated
`independent_offsets_shift3x_before_native_*` name and moves the completed native
experiment into the original path. The staging path becomes an alias so saved
provenance paths remain valid. Monitor `submitted_jobs.tsv` and `logs/` in the
native experiment. Replacement was verified complete on 2026-09-20.

The runner's optional sixth argument is the replacement target. Exact submission:

```bash
bash sim_paper/code/applications_albis/cross_modality_alignment/run_native_offsets.sh submit \
  /dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/cross_modality_alignment/native_offsets_3075_seed12345 \
  3075 12345 2025 \
  /dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/cross_modality_alignment/independent_offsets_shift3x
```

## Original independent-offset run

Run `bash submit_independent_offsets.sh` once to submit three generation/QC tasks, then STAIR and plotting after all three succeed. Jobs and logs are recorded under `data/figure_4/cross_modality_alignment/independent_offsets/`.

These datasets are exclusively for this cross-modality experiment. Figure 2 data and other analyses are neither overwritten nor redirected. Historical synchronized cross-modality results remain in the parent directory.

Generation retains the smaller Figure 2 sphere (radius 2050 µm), simulation seed 2025, weak domain mixture, modality-specific count parameters and batch sigmas (cell 1.5, bin16um 0.7, spot 0.3). Cell uses 24,207 simulated cells; bin and spot use 600,000. The native platform capture windows are retained. Per-modality rotation/translation draws are enabled by omitting `--sync-unaligned-seed`. Native modality-specific RNG offsets produce different perturbations while preserving the tissue/expression seed. Max shift is 1025 µm and max rotation 270 degrees, as in the Figure 2 generator. All modalities receive the existing QC filter.

`validate_independent_offsets.py` checks every section's fitted rotation and translation differ across modalities, verifies rigid-transform residuals, and saves `perturbation_validation.json`. Failure blocks STAIR. The dedicated STAIR job reads only these new inputs, aligning section 5 across bin16um, spot, and cell with existing STAIR settings. Its output is under `independent_offsets/STAIR/cross_tech/slice_5`; refreshed plots are under `independent_offsets/plots`.

STAIR job 35812458 completed successfully. The existing STAIR coordinate rescaling
and legacy metrics are retained in this run. Legacy per-technology independently
fitted residuals must not be interpreted as relative cross-modality alignment accuracy.

## Figure 4C reference-corrected evaluation (2026-09-20)

`reference_metrics.py` explicitly uses **bin16um** as the reference, matching the
first technology in STAIR's alignment order (spot to bin16um, then cell to spot).
For each stage (unaligned, initial, fine), fit rotation and translation using
only bin16um's observation-matched coordinates and `spatial_true`. Apply that
single transform to all three modalities. No scaling, reflection, or further
per-target fitting is allowed. Each observation is compared to its own simulated
truth; observations across modalities do not need to correspond.

The primary errors are spot and cell Euclidean RMSE in µm. The combined score is
`sqrt((RMSE_spot**2 + RMSE_cell**2) / 2)`: equal modality weight, excluding the
reference. This is not the arithmetic mean of RMSEs. The JSON also records pooled
target RMSE, per-modality medians and 95th percentiles, observation counts, and
the fitted transforms. Bin16um's residual is a reference-fit diagnostic, not a
third test modality. This is one section of one simulation seed; point-level
errors are not independent simulation replicates and no error bars are inferred.

Scores use the saved STAIR truth frame. Preprocessing scaled spot by approximately
1.0020 and cell by 1.0091 to bin16um's observed radius; it was not exactly a no-op.
The evaluation introduces no additional scale correction. These scores should
therefore be described as errors in the preprocessed coordinate frame.

Run from the project root:

```bash
MPLCONFIGDIR=/tmp/figure4c-mpl sim_paper/env/albis-tutorial/bin/python \
  sim_paper/code/applications_albis/cross_modality_alignment/reference_metrics.py \
  --input sim_paper/data/figure_4/cross_modality_alignment/independent_offsets/STAIR/cross_tech/slice_5/adata_results/Sim_CrossTech_STAIR_slice_5.h5ad \
  --outdir sim_paper/data/figure_4/cross_modality_alignment/independent_offsets/plots \
  --reference bin16um --slice 5
```

Outputs are `reference_metrics_bin16um_slice_5.{json,csv}` and
`figure4c_reference_bin16um_slice_5_{rmse,domains}.{png,pdf}`. The domain plot
separates modalities into rows and uses identical spatial limits for all panels.
Both unaligned and aligned columns use the reference-only correction. Historical
`figure4e_*` overlays use a pooled display fit and their legacy RMSE plot is
superseded by these Figure 4C files. The dedicated plotting wrapper now runs
this evaluation instead of recreating the legacy RMSE plot.

Validation: `test_reference_metrics.py` tests common rigid-motion invariance,
known target displacement, density-independent modality weighting, reference
exclusion/selection, rejection of reflection/scale correction, and invalid inputs.

## Stronger translations (2026-09-20)

The parent `cross_modality_alignment/plots/slice5_sphere_3d_domain_true.png`
belongs to the historical synchronized-offset experiment. The current independent
experiment's same-named plot is under `independent_offsets/plots/`; it now uses
bin16um-only correction for both input and result, and common spatial limits.
The vertical separation is an artificial modality display offset, not true depth.

`prepare_stronger_shifts.py` creates a separate `independent_offsets_shift3x/`
experiment from section 5 of the current QC inputs. It multiplies each modality's
original fitted XY translation by three, retaining its original rotation,
expression, observation membership and truth coordinates. It checks the resulting
transform is rigid and saves `perturbation_provenance.json`; it refuses to
overwrite an existing experiment. Original perturbation metadata is retained as
source provenance, with the changed transform recorded under `uns['stronger_shifts']`.
This experiment contains section 5 only.

`run_stronger_shifts_stair.sh` reruns STAIR and the complete plotting/metric workflow
on those changed inputs. The source experiment remains available for comparison.
The new result must be evaluated independently; source-run metrics do not apply.

Job **35816305 completed successfully** (4 min 42 s). New reference-corrected
RMSE (unaligned → STAIR fine): spot 5535.8 → 276.8 µm, cell 6036.1 → 259.8 µm;
equal-weight combined 5791.4 → 268.4 µm. The regenerated stacked comparison is
`independent_offsets_shift3x/plots/slice5_sphere_3d_domain_true.png`.

**Superseded 2026-09-20 (see "Completed replacement" above).** Job 35816305 used the
post-hoc `prepare_stronger_shifts.py` coordinate rescale, not a native ALBIS
regeneration. That data has since been archived to
`independent_offsets_shift3x_before_native_20260920T162259708032Z/` and replaced at
this same `independent_offsets_shift3x/` path by the native regeneration (jobs
35816469–471). The RMSE numbers above no longer describe what is at that path.
Current canonical numbers, from
`independent_offsets_shift3x/plots/reference_metrics_bin16um_slice_5.json`:
**spot 85.4 µm, cell 604.7 µm** (unaligned → STAIR fine). Re-read that JSON directly
before citing any number from this section — do not use the figures above.

The stacked plot now matches the alignment figure's 12 × 4.3 inch canvas,
equal column centers, 12° elevation, title row and bottom legend ("Domain True").
Each qualitative panel is centered and zoomed to fit its column; apparent sizes
across these 3D panels are therefore not a common physical display scale.
Coordinates and reference-corrected metrics are unchanged. The separate metric
domain grid retains identical spatial limits and scale across panels.
