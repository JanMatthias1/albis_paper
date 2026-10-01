# Applications of ALBIS (Figure 4)

Three downstream tasks on ALBIS 3D z-stacks (16 µm bins, spots, cells), each
using the strong-domain-mix tissue. Outputs go to `data/figure_4/<folder>/`.

## 4A: 3D z-stack alignment (`alignment/`)

A 10-slice sphere z-stack (r = 4600 µm, 1M cells) has each slice given a
random rigid offset. STAIR then aligns the slices back, using the known slice
order. Error is the RMSE of the nine non-reference sections after a rigid fit
to section 0 (no scaling or reflection).

- `generate_fig4a_tissue.sh`: simulates the tissue (prints the job id; then `bash run_stair_3D.sh <job id>`)
- `run_stair_3D.sh` → `3D_stair.py`: STAIR per modality (GPU)
- `reference_alignment_error.py`: alignment error
- `plot_figure4a.py`: before/after overlays and the error panel

## 4B: 3D spatial-domain detection (`spatial_clustering/`)

STAGATE with a 3D graph (neighbours within and across slices) vs a 2D-only
graph, scored by domain ARI, on cell, bin16um and spot × batch σ {0, 0.05,
tuned} × 3 simulation seeds (27 runs).

4B generates no data of its own: it reads Figure 3's strong-mix batch-σ data
(`data/figure_3/strong_domain_mix/batch_sigma_slide/`), so run Figure 3's
generation first ([`code/clustering`](../clustering/)):
`strong_domain_mix/generate/generate_strong_mix_{cell,bin16um,spot}.sh` (seed 2025)
and `strong_domain_mix/batch_sigma_slide/batch_sigma_slide_seeds.sh` (seeds 101, 202).

- `stagate_inputs.py`: input list (paths into the Figure 3 data)
- `run_stagate_3D.sh` → `3D_stagate.py`: training and clustering (GPU)
- `run_stagate_plots_3D.sh`, `plot_stagate_*.py`: figures

## 4C: cross-modality alignment (`cross_modality_alignment/`)

One ALBIS call builds a single tissue that is captured as bin16um, spot and
cell. Each modality gets its own rigid offset (up to 1.5 × the sphere radius).
STAIR aligns spot and cell to bin16um. Error is the RMSE against the bin16um
reference, plus how similar the genes are across modalities after alignment.

- `strong_domain_mix/run_strongmix_offsets.sh submit`: generation → STAIR → metrics → plots
- `cross_tech_stair.py`, `reference_metrics.py`, `evaluation_gene_metrics.py`: alignment and metrics (shared)
- `plot/`: figures

### Cropped bin/spot variant

The same tissue and offsets, but with bin and spot cropped to a 2221 µm square
capture area (Visium's 6.5 mm scaled to the smaller sphere) while cell stays
uncropped. Alignment uses `cross_tech_stair.py --no-rescale`, so every modality
keeps its true size and the bin/spot squares sit inside the cell section.

- `run_cropped_no_rescale.sh submit`: generation → STAIR → metrics → plots

## Simulation parameters

All three panels use the strong domain mix: 6 spatial domains, 8 cell types,
each domain enriched for two cell types (0.30 each, 0.067 for the rest);
556 genes (80 markers per type, 25% shared markers, 10% noise genes);
10 slices along Z; cell radius lognormal (mean 7.5 µm, σ 0.28, limited to
4–14 µm), the same in every panel; no overlapping cells.

| | 4A | 4B (Figure 3 data) | 4C |
|---|---|---|---|
| Tissue | one shared tissue for all modalities | a separate tissue per modality | one shared tissue for all modalities |
| Sphere radius | 4600 µm | 2050 µm | 2050 µm |
| Cells | 1,000,000 | cell 24,207; bin/spot 600,000 | 600,000 |
| Gene mean (lognormal μ, σ) | −2.5, 0.7 | cell −2.5, 0.7; bin16um −2.5, 0.7; spot −2.25, 1.0 | −2.5, 0.7 |
| Dispersion θ (jitter) | 2.0 (0.6) | cell 0.40 (0.15); bin16um 2.0 (0.6); spot 0.25 (0.10) | 2.0 (0.6) |
| Batch σ | 0.7 | {0, 0.05, tuned}: cell 1.5, bin16um 0.7, spot 0.3 | 0.7 |
| Per-slice offset (max shift, max rotation) | 6900 µm (1.5 R), 270° | 1025 µm (0.5 R), 270° (not used: STAGATE uses the aligned coordinates) | 3075 µm (1.5 R), 270° |
| Domain-boundary fuzz | 230 µm | 102.5 µm | 102.5 µm |
| Simulation seed | 2025 | 2025, 101, 202 | 2025 |

4A and 4C use the Figure 3 bin16um expression settings for all three
modalities, so spot and cell there are not calibrated to their Figure 2 real
references. Capture: bins 16 µm and spots 55 µm diameter at 100 µm spacing,
both within a 6.5 mm square (Visium / Visium HD); cells within a 12 × 24 mm
window (Xenium). The cropped 4C variant instead captures bin and spot within a
2221 µm square. 4A's cell density is about 1/7 of Figures 2 and 3 (illustrative).
Each tissue's full settings are saved as `config.json` next to its data.

## Environments

STAIR (`STAIR-tools`, GPU, with R/mclust) for 4A and 4C;
`env/stagate-pyg` (GPU) for 4B.
