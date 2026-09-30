# ALBIS 16um-parameter tissue: spot RCTD and cell-type clustering

Three newly generated simulations (seeds 2025, 101, 202), including a fresh
seed2025. Exact Figure 5A / Figure 4 strong_domain_mix_shift3x config; only the
tissue seed varies. All modalities generated together natively. 600k cells,
556 genes, 8 types, 6 domains, 10 slices, theta2, theta_jitter0.6,
base_gene_lognormal[-2.5,0.7], batch_sigma0.7. Spatial perturbation seed12345
remains fixed. 16um is the bin-derived parameter set; query spots are circular,
diameter55um, spacing100um.

Each seed directory contains config.json, data/{cell,bin16um,spot}/,
generation_status.json, inputs/ (native counts, excluding empty/zero/unassigned),
results/ (RCTD predictions, object, metrics, scatter), rctd.log, and clustering/.
Reference uses matching cells from the same tissue and seed; not independent.
RCTD uses the unchanged shared weak_domain/run_rctd_spot.R with seed0,
full mode and default QC; truth is native captured-molecule fractions.
All reference cells are supplied; RCTD retains its default reference handling.

Code:
- generate.py --seed SEED: native full shared-tissue simulation from saved config.
- prepare_inputs.py --seed SEED: validate count scale and prepare RCTD inputs.
- run.sh: generation + input preparation + RCTD, array0-2.
- cluster_cells.py --seed SEED / run_clustering.sh: no-Harmony PCA30, Gaussian
  K256, Leiden targeting8, fixed clustering seed0. Positive-total cells only;
  raw generated cells retain batch_sigma0.7. Waits for the current clustering
  experiment and this generation/RCTD array. No old results are replaced.
- plot_deconvolution.py: single-axis three-seed comparison to existing SPIDER
  556/2000-gene runs, PNG/PDF/SVG, per-seed/mean-SD CSVs and provenance in figures/.
- summarize.py: ALBIS seed metrics and Slice5 true/estimated spatial panels.
- finish.sh: validate/summarize and plot after RCTD succeeds for all3 seeds.

jobs.json in the data directory records submitted job IDs and dependencies.
Scripts run in project environments: sim_paper/env/albis-tutorial (generation,
clustering), comparison_methods/env/analysis (preparation/plots), sim_paper/env/rctd.
Existing output directories are protected against accidental reruns.
Differences from old Figure3 deconvolution: old cells used batch_sigma1.5,
spots0.3 and a fixed query; these new replicates use a shared tissue per seed,
batch_sigma0.7, and both query and reference vary.
