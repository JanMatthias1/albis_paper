# Applications of ALBIS (Figure 4)

Three downstream tasks on ALBIS 3D z-stacks (16 µm bins, spots, cells), each
using the strong-domain-mix tissue. Outputs go to `data/figure_4/<folder>/`.

## 4A: 3D z-stack alignment (`alignment/`)

A 10-slice sphere z-stack (r = 4600 µm, 1M cells) has each slice given a
random rigid offset. STAIR then aligns the slices back, using the known slice
order. Error is the RMSE of the nine non-reference sections after a rigid fit
to section 0 (no scaling or reflection).

- `run_stair_3D.sh` → `3D_stair.py`: STAIR per modality (GPU)
- `reference_alignment_error.py`: alignment error
- `plot_figure4a.py`: before/after overlays and the error panel

## 4B: 3D spatial-domain detection (`spatial_clustering/`)

STAGATE with a 3D graph (neighbours within and across slices) vs a 2D-only
graph, scored by domain ARI. The data are the Figure 3 strong-mix batch-σ
inputs: cell, bin16um and spot × batch σ {0, 0.05, tuned} × 3 simulation
seeds (27 runs).

- `stagate_inputs.py`: input list
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

## Environments

STAIR (`STAIR-tools`, GPU, with R/mclust) for 4A and 4C;
`env/stagate-pyg` (GPU) for 4B.
