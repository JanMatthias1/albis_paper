# Current cell-type comparison

600,000 cells across 10 slices per realization; simulation seeds 2025, 101, 202.
ALBIS: 556 genes. scCube and SPIDER: 556 and 2,000 genes. Total: 15 runs.

ALBIS uses the uncropped Figure 5A / Figure 4 strong_domain_mix_shift3x configuration (16 µm reference). Only batch_sigma=1.5, simulation seed, and cell-only output selection change. SPIDER uses native simulate_10X_3d with eight equally abundant cell types, self-neighbor probability 0.7 and the Figure 5A iteration settings; this is separate from gyrus domain recovery. scCube uses grid size 8, delta 2, lambda 0.75 and 200 VAE epochs. Each gene-count/seed has a freshly generated 10,000-cell Splatter reference shared by scCube and SPIDER; each VAE is trained anew.

`cluster.py` directly imports `../../../clustering/step01_pca_harmony.py` and `step02_leiden_resolution_sweep.py` via their absolute project paths. All methods use positive-total cells, normalize_total(10000), log1p, scale(max_value=10), 30 PCs, Harmony(slice_id), 15 neighbors and Leiden with fixed seed 0. Resolution targets eight clusters, never maximum ARI. The actual achieved count is reported.

## Code and execution

- `prepare_experiment.py`: freeze settings and manifests in the new experiment directory; refuses overwrite.
- `run_task.py` / `run_array.sh`: generate six reference pools, then generate and cluster 15 datasets using the appropriate environments.
- `generate_albis.py`: native ALBIS generation.
- `cluster.py`: common analysis, per-cell labels, embeddings and ARI.
- `aggregate.py` / `aggregate.sh`: status and metric tables; final boxplot only when all 15 metrics exist.

Competitor generation reuses `../../overview/generate.py` with `cell_only=true`. Analysis uses `sim_paper/env/albis-tutorial`; competitor and Splatter environments remain under `comparison_methods/env`.

Submit pool array 0-5, then run array 0-14 with afterok dependency on pools, then aggregation with afterany dependency on runs. Job identifiers are recorded in the experiment's jobs.json. Review failures in per-run logs before retrying; generation and clustering refuse existing output directories to prevent accidental overwrites.

## Outputs

`sim_paper/data/figure_5/clustering/cell_type/experiment_600k/`:

- `protocol.json`, `pools.json`, `runs.json`, `jobs.json`: design and execution records.
- `pools/genes*_seed*/`: contracts, Splatter reference matrices and generation logs.
- `runs/{method}_genes{n}_seed{seed}/`: settings, generation/clustering logs, `{method}/cell.h5ad` and manifest.
- Each run's `clustering/`: `metrics.json`, `labels.csv.gz`, `embeddings_labels.h5ad` (embeddings/labels/graph, no duplicated expression).
- `status.csv`, `metrics.csv`, `cell_type_ari_boxplot.{png,pdf,svg}` after aggregation.

ALBIS includes an imposed batch effect; competitors do not have an equivalent imposed batch effect. Gene-count conditions regenerate expression rather than subset one matrix. Three seed points measure simulation variation. These results compare recoverability under the stated settings, not an identical-noise benchmark.

## Archive

Previous code and results are preserved in the respective code/data `archive/pre_600k_20260927/` directories, with `moved_files.json` inventories. Historical scripts retain their original contents and may contain pre-archive paths; they are not the active entry points.

## SPIDER graph repair pilot

`repair_graph.py` groups exact raw-expression profiles by slice (SHA256 of canonical CSR entries plus slice_id), averages their saved Harmony embeddings, constructs a 15-neighbor graph on unique profile/slice groups, runs the existing known-eight Leiden search, and broadcasts labels to all original cells. Grouping never uses cell-type truth. Results are separate under `clustering_duplicate_aware/`; the original analysis is retained. This changes profile weighting in graph construction and needs validation/common application before replacing the main comparison. Pilot job 35989515 uses SPIDER 556 genes, seed 2025.

Validated follow-up: `--group-across-slices` collapses identical raw profiles across slices, averages their saved Harmony embeddings, builds the 15NN graph, and broadcasts labels. In the SPIDER 556/seed2025 pilot this yielded 10,000 unique profiles, 8 connected components, 8 Leiden clusters and ARI=1.0 (the within-slice-only pilot still yielded 10,000 clusters). No labels were used to construct groups or the graph. This is unweighted unique-profile graph clustering, with ARI evaluated over all original cells; it changes graph weighting and must be reported. No simulation/noise was changed.

Results go to `clustering_unique_profiles/` within each run. If no duplicate profiles exist, the original analysis is retained. Common-policy array 35990335 runs after original array 35986639; summary job 35990340 writes to experiment_600k/clustering_unique_profiles/. Original metrics are preserved, and summaries refuse to present count-mismatched runs as a completed comparison.

## Current direction: full-cell KNN, no duplicate collapsing

User rejected the duplicate-collapse approach. Cancelled queued jobs 35990335 and 35990340; pilot outputs remain diagnostics only. New code `knn_graph.py` / `knn_graph.sh` retains all original cells and existing PCA/Harmony embeddings. Pilots 35990644 test UMAP K=128, UMAP K=256, Gaussian K=256 on SPIDER 556/seed2025. Selection job 35990671 uses connectivity only (first candidate with <=8 components), never ARI. Gaussian changes edge weighting as well as K; the chosen method is recorded in knn_protocol.json. If no candidate passes, selection fails explicitly.

The proposed 15-run follow-up submission was rejected by the user and was NOT submitted. Review pilot diagnostics before any further full-array submission. Deconvolution job 35989955 is unaffected.

## Paired Harmony diagnostic (before deciding graph changes)

User requested causal diagnosis before changing the pipeline. Automatic K selection job 35990671 was cancelled. Diagnostic array 35990742 runs `diagnose_harmony_graph.py` for SPIDER 556 and 2000 genes, seed 2025. It rebuilds K=15 UMAP-weighted Euclidean graphs from saved pre-Harmony and post-Harmony PCA, keeping all 600,000 cells/order and random seed 0 unchanged. It reports component sizes, degree distributions and whether components span slices, alongside the original saved post-Harmony graph. Outputs: experiment_600k/harmony_graph_diagnostic/. No graph settings or clustering results are adopted by this diagnostic.

## Requested Harmony policy update

scCube/SPIDER now use pre-Harmony PCA; ALBIS retains Harmony(slice_id). `cluster.py` skips Harmony for competitors and writes separate `clustering_no_harmony/` outputs. Existing all-Harmony outputs remain historical and are not relabelled. Default aggregation now requires the new competitor outputs and will not silently mix in historical Harmony results. New full-KNN pilots use pre-Harmony PCA for competitors and separate output/protocol paths, preventing reuse of the previous post-Harmony selection. No new graph policy has been adopted or reruns submitted with this change; K=15 remains the baseline until a graph is chosen.

Completed earlier post-Harmony SPIDER pilot connectivity: UMAP128 8,342 components; UMAP256 6,739; Gaussian256 8 components of 75,000 cells each. These diagnostics favor testing Gaussian256 on the newly requested pre-Harmony basis; they do not establish pre-Harmony performance or ARI.

## Active trial: all methods without Harmony, batch-zero ALBIS

User-selected ALBIS source is Figure 3 `strong_domain_mix/batch_sigma_slide/cell/bs0{,_seed101,_seed202}` QC files (about 24k cells, 556 genes). `batchzero_trial.py` and `.sh` run a separate 15-condition trial under `data/figure_5/clustering/cell_type/batchzero_noharmony_gauss256/`: no Harmony for any method, K=256 Gaussian graph for all, no profile collapsing, fixed seed 0 and known-eight Leiden search. Competitors retain the 600k generated datasets and saved pre-Harmony PCA. Job array 35992151; dependent summary 35992155. The cell-count mismatch is recorded. This supersedes earlier primary-pipeline plans; historical scripts/results are retained for provenance.
