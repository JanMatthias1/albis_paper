# ALBIS: simulation of multi-resolution and multi-dimensional spatial transcriptomics data

[![PyPI](https://img.shields.io/pypi/v/albis)](https://pypi.org/project/albis/)
[![ALBIS package](https://img.shields.io/badge/GitHub-JanMatthias1%2FAlbis-blue)](https://github.com/JanMatthias1/Albis#readme)

This repository contains the analysis and figure code for the ALBIS manuscript.
The ALBIS package itself is developed at
[JanMatthias1/Albis](https://github.com/JanMatthias1/Albis#readme).

**Paper:** link to be added\
**Authors:** Jan Matthias\*, Jianing Yao\*, Stephanie C. Hicks\
\* These authors contributed equally.

## Abstract

The rapid adoption of spatial transcriptomics (ST) technologies has been accompanied by substantial growth in computational methods, creating an increasing need for systematic and reliable benchmarking, particularly as analyses increasingly extend from individual tissue sections to three-dimensional (3D) tissue reconstruction. Existing benchmarks often rely on a small and similar set of experimentally derived datasets because they offer reliable ground truth labels, which are hard to find in practice. Simulated or synthetic data provide an alternative, yet most existing simulators either adapt single-cell RNA-sequencing models to spatial coordinates or generate individual sections in a two-dimensional (2D) space at predefined resolutions. More recent approaches can simulate 2D and 3D tissues, but generally do not derive multiple spatial resolutions from the same molecular realization.

Here, we introduce ALBIS, a statistically grounded and highly customizable framework that generates synthetic ST data from an *in silico* 3D tissue. ALBIS defines spatially structured cell populations with corresponding markers for gene expression, then initializes individual mRNA molecules in 3D space to create cell-level information with a ground truth 3D coordinate for each molecule and cell. From this common realization, ALBIS cuts the 3D tissue into multiple 2D tissue sections with cell-, bin-, and spot-level measurements, enabling direct comparison across spatial resolutions. Properties of the simulated data, including spatial domains, cell-type composition, gene-expression effects, batch effects, transcript spillover, measurement resolution, and tissue scale, can be systematically varied, enabling rigorous stress-testing of computational methods under known and increasingly challenging conditions.

## Installing ALBIS

```bash
pip install albis
```

## Figures

Workflows and run commands are described in the READMEs inside each folder.

| Figure | Content | Code |
| --- | --- | --- |
| 1 | The 3D tissue sphere and how it is captured as cells, bins and spots | [`code/sphere_figure`](code/sphere_figure/) |
| 2 | Count distributions of simulated vs real data (Xenium, Visium HD, Visium) | [`code/count_distribution`](code/count_distribution/) |
| 3 | Spatial-domain and cell-type recovery, batch effects, spot deconvolution | [`code/clustering`](code/clustering/) |
| 4 | Applications: 3D alignment, 3D spatial domains, cross-modality alignment | [`code/applications_albis`](code/applications_albis/) |
| 5 | Comparison with other spatial simulators (scCube, SPIDER) | [`code/comparison_methods`](code/comparison_methods/) |

Shared code:

| Path | Content |
| --- | --- |
| [`code/data`](code/data/) | generates the ALBIS simulations used across figures |
| [`code/real_data_qc`](code/real_data_qc/) | QC of the real 10x reference datasets |
| [`code/manuscript_style.py`](code/manuscript_style.py) | shared plot style (fonts, colors) |
| [`env`](env/) | scripts that create the conda environments (ALBIS, BANKSY, STAGATE, RCTD) |

## Running the code

- Simulated and real data are not included in this repository. Scripts read and
  write under `data/`, and the real datasets are public 10x Genomics downloads.
- Most workflows run as SLURM jobs, and several need a GPU (STAIR, STAGATE).
  Paths and SLURM settings come from the original computing cluster, so adapt
  them before running.
- The main environment is created with `env/create_tutorial_env.sh`; tools with
  conflicting dependencies (BANKSY, STAGATE, RCTD) each have their own.

## Contact

Questions or comments about the manuscript:
Stephanie C. Hicks (corresponding author), [shicks19@jhu.edu](mailto:shicks19@jhu.edu)
