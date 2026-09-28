# Figure finalization audit — 2026-09-19

## Update — 2026-09-26
- **Figure 4C source changed** to `data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x_cropped_bin_spot/`
  (strong mix, shared tissue, bin16um + spot cropped to a 2221 µm square).
  Spot RMSE 5041.5 → 16.7 µm, cell 4441.3 → 18.4 µm; gene PCC 0.183 (STAIR) vs
  0.181 (truth). Replaces `independent_offsets_shift3x/`. Caption caveats (cropping
  eases alignment; gene-PCC ceiling depends on the cropped overlap; single section/seed)
  are listed in `code/applications_albis/cross_modality_alignment/README.md`.
- **STAGATE (4B) is being rerun**: the 2026-09-18/19 results used strong-mix data
  regenerated on 2026-09-23. Jobs 35933261–87; 26/27 done, bin16um bs0.7 seed202
  waiting for an H100. The "27/27 plotting" line below is from the old data; replot
  once the last job finishes.
- **Figure 3 summary** (`banksy_vs_pca_recovery.png`) now shows mean ± SD over
  simulation seeds 2025/101/202 for every bar (BANKSY λ fixed per modality/target).

## Update — 2026-09-20
- STAIR 35812458 completed successfully; no resubmission needed.
- Figure 4C now has a separate reference-corrected evaluator, explicitly anchored
  to bin16um. Spot RMSE: 3296.6 → 105.7 µm; cell: 3130.5 → 325.9 µm.
- Six analytic synthetic tests passed. New Figure 4C plots and JSON/CSV use
  `figure4c_reference_bin16um_slice_5_*` / `reference_metrics_bin16um_slice_5.*`
  in the independent-offset experiment's `plots/` directory. Legacy metrics remain
  for provenance, not final scoring. Panel selection/assembly remains open.
- Scores use the saved preprocessed truth frame: spot/cell were scaled by about
  1.0020/1.0091 before alignment. See the dedicated experiment README for methods.

The audit below records the prior day's state.

Reviewed plotting style guide; count-distribution, clustering, real-QC, STAGATE,
alignment and independent-offset workflow READMEs; checked current output tables,
completion logs and cross-modality queue status. Historical `figure.md` and
`DATA_VERSIONS.md` contain superseded panel labels/status and are not final manifests.

## Verified
- Figure 1 palette refresh jobs 35812471–35812473 finished successfully.
- STAGATE plotting array finished 27/27 runs. Updated batch summary validates
  all three seeds for all nine conditions (18 method/condition summaries).
- New cross-modality generation finished for bin16um, spot and cell. Retained
  observations: 365785, 10034, 24207. Independent rotations/translations verified
  across modalities in all ten sections; provenance JSON saved in experiment root.
- STAIR job 35812458 is PENDING (Priority) at audit; no finished new alignment yet.
- Real-QC supplementary figure contains nine panels: seven tissues/platform
  samples, with two HD tissues each plotted at two resolutions. Both HD resolutions
  are not independent samples. QC membership checked against saved filtered files.
- Figure 4A reference RMSE panel uses cell_r6000, bin16um and spot; reference 0
  excluded from scored distributions. Detailed methods file is present.

## Remaining before final export
1. Finish independent-offset STAIR, review outputs, add/validate reference-corrected
   cross-modality metric. Use only the dedicated experiment; parent-directory
   synchronized-offset plots are historical. Select the reference modality explicitly.
2. Lock the actual selected panel files and final numbering (4A alignment,
   4B STAGATE, 4C cross-modality). File names still include old 4C/4E labels.
   Record all choices in a panel manifest; do not infer inclusion of historical SVG panels.
3. Replot selected files with latest code where updates were applied only to
   examples. In particular, STAGATE larger spheres/elevation were rerendered for
   spot_strongmix_bs0_seed101 only; Figure 3 layout refinements also affected subsets.
4. Captions: simulation seeds versus training seeds; mean ± sample SD; RCTD
   reference averaging and violin display cap; 4A larger cell geometry and
   reference-only metric versus whole-stack fit used for qualitative display;
   nine QC panels versus seven distinct sample/platform datasets.
5. Assemble final main/supplementary figures and inspect at final publication size.
   Check fonts, whitespace, legend identity, chosen shift variant in Figure 1,
   raster/export quality and full script/input/output provenance.

No broad retraining of Figures 1–3 or STAGATE is indicated by this audit.
The current scientific computation still pending is the dedicated cross-modality run.
