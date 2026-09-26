# Figure 4C: cross-modality alignment

STAIR aligns bin16um, spot and cell captures of the same section (slice 5), each
given its own random rigid offset. Outputs go to
`data/figure_4/cross_modality_alignment/<experiment>/`.

## Current Figure 4C source (chosen 2026-09-26)

`data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x_cropped_bin_spot/`:
strong domain mix, one shared tissue, shift3x offsets, with bin16um **and** spot
cropped to a 2221 µm square capture area (Visium's 6.5 mm scaled by the sphere
shrink 2050/6000); cell uncropped. Panels: `plots/figure4c_reference_bin16um_slice_5_{domains,rmse}.png`,
`plots/slice5_sphere_3d_domain_true.png`, `plots/figure4e_overlay_2d_*`,
`gene_similarity/gene_similarity_ridges_average_*.png`.

Section 5, bin16um reference: spot RMSE 5041.5 → 16.7 µm, cell 4441.3 → 18.4 µm;
gene PCC median 0.183 after STAIR vs 0.181 at ground truth.

Caption notes:
- Cropping makes alignment much easier than on the round, uncropped tissue
  (`strong_domain_mix_shift3x/`: spot 353 µm, cell 394 µm, ~20-30° leftover
  rotation); the square capture corners probably pin the rotation (not tested).
- Gene PCC is computed where bin16um overlaps, i.e. the central square (mostly
  domain D5), so its ceiling is lower (truth 0.18) than uncropped (0.48); compare
  STAIR with truth within a run, not across runs.
- One section, one simulation seed: illustrative, not a robustness claim.

Replaced as the headline: `independent_offsets_shift3x/` (weak mix, one ALBIS
call per modality). It and the other runs stay on disk.

| Path | Contents |
|---|---|
| `strong_domain_mix/` | Strong-domain-mix, shared-tissue experiment: one ALBIS call builds a single tissue for all three modalities. `bash run_strongmix_offsets.sh submit [OUTDIR]` -> `strong_domain_mix_shift3x/`. Add `--crop-modalities bin[,spot] --crop-window-um 2221` to crop bin/spot to a square capture area on the same tissue -> `strong_domain_mix_shift3x_cropped_bin/`, `..._cropped_bin_spot/` |
| `weak_domain_mix/` | Canonical (weak-mix) native-offset experiment: one ALBIS call per modality. `bash run_native_offsets.sh submit ...` -> `independent_offsets_shift3x/`. Details: `weak_domain_mix/README_independent_offsets.md` |
| `plot/` | Plotting. `plot_independent_offsets.py --base OUT` is the driver both pipelines call; it imports the other three plot modules |
| `misc/` | Retired code, kept for provenance only. Paths inside these files were not updated after the 2026-09-25 reorganization, so they will not run as-is |

Shared scripts (top level, called by both pipelines):

- `cross_tech_stair.py`: STAIR cross-technology alignment (runs in the STAIR conda env).
- `reference_metrics.py`: bin16um-reference rigid correction and held-out RMSE, the Fig4C RMSE numbers. `test_reference_metrics.py` holds its unit tests.
- `evaluation_gene_metrics.py`: all-gene cross-modality grid similarity, the headline Fig4C metric. `run_evaluation_gene_metrics.sh` runs it standalone.

`metrics.json` from `cross_tech_stair.py` includes `stair_fine_rmse_per_technology_um`.
That field is always about 0 because it refits each modality on its own. Don't cite it.
