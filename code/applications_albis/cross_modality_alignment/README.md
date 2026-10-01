# Figure 4C: cross-modality alignment

STAIR aligns bin16um, spot and cell captures of the same section (slice 5), each
given its own random rigid offset. Outputs go to
`data/figure_4/cross_modality_alignment/<experiment>/`.

## Current Figure 4C source (chosen 2026-10-01, albis 0.1.2)

`data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x/`: strong domain
mix, one shared tissue (r = 2050 um, 600k cells, bin16um expression settings for all
three modalities), shift3x offsets (max shift 3075 um = 1.5 R, 270 deg), **uncropped**.
`bash strong_domain_mix/run_strongmix_offsets.sh submit` (defaults). The same folder
is Figure 5A's reproducibility reference. Expect larger residual error than the
cropped run below (on albis 0.1.1: spot 353 um, cell 394 um, ~20-30 deg rotation).

Figure 4A has its own larger tissue from the same generator
(`--sphere-r-um 4600 --n-cells 1000000` -> `strong_domain_mix_shift3x_r4600/`);
see `../alignment/README.md`.

## Previous Figure 4C source (2026-09-26, superseded; not rerun on albis 0.1.2)

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
| `plot/` | Plotting. `plot_independent_offsets.py --base OUT` is the driver the pipeline calls; it imports the other three plot modules |

Retired code lives outside this folder, in `code/misc/applications_albis/cross_modality_alignment/`:
`weak_domain_mix/` (the replaced weak-mix native-offset experiment, one ALBIS call per
modality -> `independent_offsets_shift3x/`; moved 2026-10-01, not rerun on albis 0.1.2) and
`misc/` (older helpers; paths not updated after the 2026-09-25 reorganization, so they will
not run as-is). All Figure 4 data from albis <= 0.1.1 is in `data/figure_4/archive_albis0.1.1_20261001/`.

Shared scripts (top level, called by both pipelines):

- `cross_tech_stair.py`: STAIR cross-technology alignment (runs in the STAIR conda env).
- `reference_metrics.py`: bin16um-reference rigid correction and held-out RMSE, the Fig4C RMSE numbers. `test_reference_metrics.py` holds its unit tests.
- `evaluation_gene_metrics.py`: all-gene cross-modality grid similarity, the headline Fig4C metric. `run_evaluation_gene_metrics.sh` runs it standalone.

`metrics.json` from `cross_tech_stair.py` includes `stair_fine_rmse_per_technology_um`.
That field is always about 0 because it refits each modality on its own. Don't cite it.
