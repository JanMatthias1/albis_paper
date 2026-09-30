# Preliminary cell-type clustering

SPIDER and scCube use the completed native 10,000-cell, 2,000-gene datasets in
`data/comparison_methods/overview_full_35831924`. All cells are validated for
finite nonnegative expression, positive totals, unique IDs, eight known types,
and valid slice labels. Input data are not changed; scCube continuous values
are used directly, without the rounding required by the RCTD pilot.

This wrapper reuses the same implementation as
`code/clustering/strong_domain_mix/pca_harmony/strongmix_celltype_panel.sh`:

1. All genes, normalize total to 10,000, log1p, scale with clipping at 10.
2. PCA with 30 components, Harmony using slice_id, 15-neighbor graph.
3. Leiden resolution search over 0.01–3.0, at most 15 bisection trials,
   targeting the eight known cell types. Select the closest cluster count,
   not the highest ARI. Report achieved count if eight cannot be reached.
4. ARI against cell_type_true, before/after Harmony diagnostics, final UMAP
   and true-versus-predicted contingency plots. Random state is the original
   scripts' default 0.

The PCA and plotting scripts are invoked directly. The existing resolution
script is imported with GROUND_TRUTH_COLS restricted to cell_type_true, because
competitors do not have Albis domain labels. The search function is unchanged;
no artificial domain_true column is added. Each run records input and shared
script checksums in input_audit.json and retains every resolution trial.

## Submit

From sim_project (create this folder's logs/ before submission):

```bash
sbatch --output=sim_paper/code/comparison_methods/clustering/cell_type/logs/spider_%j.out sim_paper/code/comparison_methods/clustering/cell_type/run_pilot.sh spider
sbatch --output=sim_paper/code/comparison_methods/clustering/cell_type/logs/sccube_%j.out sim_paper/code/comparison_methods/clustering/cell_type/run_pilot.sh sccube
```

Each job requests 4 CPUs, 32 GB, and four hours. Outputs are isolated under
`data/comparison_methods/clustering/cell_type/pilot_METHOD_JOBID/`.
Read `ari_recovery/ari_summary_cell.json` for ARI and resolution trials.
Initial jobs: SPIDER 35973359; scCube 35973360.

## Interpretation

This reproduces the Albis analysis procedure, not a controlled simulator
comparison. Competitors have 2,000 genes versus Albis's 556; marker strength,
expression distributions, geometry and slice composition differ. There is no
explicit Albis-style slice batch effect in these competitor outputs. Harmony
may therefore also remove genuine slice-associated biology; before/after
plots help inspect that, but no PCA-only ARI control is included in this pilot.
Known-k tuning uses the true category count and should be disclosed as such;
ARI evaluates assignments and is not used to select the resolution. No broad
PCA/neighborhood parameter optimization or independent-seed replication is
included in this preliminary run.

## Native 556-gene comparison

Jobs SPIDER 35973656 and scCube 35973657 use the freshly regenerated datasets
in `data/comparison_methods/overview_native556_20260926/`, not subsets of the
2,000-gene outputs. They retain the same normalization, 30 PCs, Harmony by
slice, 15 neighbors and known-eight-type Leiden resolution tuning.
The wrapper accepts CELL_CLUSTER_SOURCE, CELL_CLUSTER_PREFIX and
CELL_CLUSTER_N_GENES; these submissions set native556 and require exactly
556 unique input features. Outputs are `native556_METHOD_JOBID/` beneath
`data/comparison_methods/clustering/cell_type/`.

They wait for jobs 35973632/35973633 to end using afterany, because those jobs
combine generation with RCTD: a later deconvolution failure need not block
clustering of valid generated cell data. Missing or invalid cell inputs cause
clustering to fail explicitly. The original 2,000-gene results remain separate.
Changing the generated gene count can also change expression distributions;
this comparison is not purely removal of features from the same expression matrix.
