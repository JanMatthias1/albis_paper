# Manuscript figures

For the next session, start with the **Next session — final figure checklist** in
[figure.md](figure.md). The detailed verified status and remaining work are in
[FIGURE_FINALIZATION_AUDIT.md](FIGURE_FINALIZATION_AUDIT.md).

Cross-modality STAIR job **35812458** completed successfully. Reference-corrected
Figure 4C evaluation now uses **bin16um** as the reference; see
[the experiment README](code/applications_albis/cross_modality_alignment/README_independent_offsets.md).
The experiment and new `figure4c_reference_bin16um_slice_5_*` plots are isolated
under `data/figure_4/cross_modality_alignment/independent_offsets/`.

Latest requested stronger-displacement variant: `independent_offsets_shift3x/`.

Native regeneration and replacement completed successfully on 2026-09-20:
generation **35816469** (all three tasks), alignment/plots/metrics
**35816470**, automatic replacement **35816471**, all exit code 0.
All ten sections were regenerated; alignment and evaluation cover section 5.
The previous version is archived as `independent_offsets_shift3x_before_native_20260920T162259708032Z/`.
Logs and job IDs: `data/figure_4/cross_modality_alignment/native_offsets_3075_seed12345/`.

Historical STAIR job **35816305** used posthoc 3× translations and is now archived.
The current `plots/slice5_sphere_3d_domain_true.png` uses the native regenerated inputs.
The parent `cross_modality_alignment/plots/` still contains historical synchronized results.
