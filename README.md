# ALBIS: simulation of multi-resolution and multi-dimensional spatial transcriptomics data

The rapid adoption of spatial transcriptomics (ST) technologies has been accompanied by substantial growth in computational methods, creating an increasing need for systematic and reliable benchmarking, particularly as analyses increasingly extend from individual tissue sections to three-dimensional (3D) tissue reconstruction. Existing benchmarks often rely on a small and similar set of experimentally derived datasets because they offer reliable ground truth labels, which are hard to find in practice. Simulated or synthetic data provide an alternative, yet most existing simulators either adapt single-cell RNA-sequencing models to spatial coordinates or generate individual sections in a two-dimensional (2D) space at predefined resolutions. More recent approaches can simulate 2D and 3D tissues, but generally do not derive multiple spatial resolutions from the same molecular realization.

Here, we introduce ALBIS, a statistically grounded and highly customizable framework that generates synthetic ST data from an *in silico* 3D tissue. ALBIS defines spatially structured cell populations with corresponding markers for gene expression, then initializes individual mRNA molecules in 3D space to create cell-level information with a ground truth 3D coordinate for each molecule and cell. From this common realization, ALBIS cuts the 3D tissue into multiple 2D tissue sections with cell-, bin-, and spot-level measurements, enabling direct comparison across spatial resolutions. Properties of the simulated data, including spatial domains, cell-type composition, gene-expression effects, batch effects, transcript spillover, measurement resolution, and tissue scale, can be systematically varied, enabling rigorous stress-testing of computational methods under known and increasingly challenging conditions.

**Authors:** 

**Paper:** Link to be added.

## The ALBIS package

ALBIS is available on PyPI as [`albis`](https://pypi.org/project/albis/):

```bash
pip install albis
```

Source code and package documentation are on GitHub:
[JanMatthias1/Albis](https://github.com/JanMatthias1/Albis#readme).
This repository contains only the manuscript analyses that use it.

## Code for the manuscript figures

This repository contains the analysis and plotting code for the ALBIS manuscript. Code is organized by analysis; some directories contribute to more than one figure.

| Figure | Code directory |
| --- | --- |
| Figure 1 | [code/sphere_figure](code/sphere_figure/) |
| Figure 2 | [code/count_distribution](code/count_distribution/) |
| Figure 3 | [code/clustering](code/clustering/) |
| Figure 4 | [code/applications_albis](code/applications_albis/) |
| Figure 5 | [code/comparison_methods](code/comparison_methods/) |

Analysis-specific README files and shell scripts provide details on individual workflows. Some scripts use paths and SLURM settings from the original computing environment; adapt these to your installation before running them.

## Contact

Questions or comments related to the manuscript:

- Stephanie C. Hicks (corresponding author): [shicks19@jhu.edu](mailto:shicks19@jhu.edu)
