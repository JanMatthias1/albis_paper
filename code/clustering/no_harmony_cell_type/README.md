# Figure 3: cell-type recovery without Harmony

Code: `sim_paper/code/clustering/no_harmony_cell_type/`.
Outputs: `sim_paper/data/figure_3/no_harmony_cell_type/`.
The existing domain experiment remains in `data/figure_3/no_harmony_domain/`.

Uses the same 18 batch-zero QC datasets (556 genes, 10 slices), strong/weak domain
mixtures, cells/16 µm bins/spots, and simulation seeds 2025, 101, 202. The cell
datasets contain 24,207 observations; this is separate from Figure 5's 600k runs.
Bins/spots are evaluated against their stored dominant cell-type label.

Sweep seven BANKSY lambdas (0,.1,.2,.3,.5,.8,1) on strong seed2025, separately
for each modality. Fixed k_geom: 60 cells, 100 bins, 8 spots. Choose maximum
cell-type ARI among runs recovering exactly eight clusters; exact ties choose
the smaller lambda. This is a target-count evaluation, not unsupervised selection
of the number of types. Seed2025 participates in tuning and reporting.

Reuse existing no-Harmony embeddings and K15 graphs when the data and BANKSY
settings match. Reclustering targets eight cell types instead of six domains.
Missing selected-lambda embeddings use the unchanged domain embedding routine.
BANKSY preprocessing: normalize total to 10,000 without log, max_m=1,
scaled_gaussian, audited within-slice spatial edges, PCA30. Expression-only:
normalize10000, log1p, scale clipped10, PCA30. Both use K15 UMAP-weighted neighbors
and igraph Leiden, seed0, n_iterations2, resolution .01–3, up to15 trials.
No Harmony is applied. Existing Figure 3 tuned-batch panels remain unchanged.

`pipeline.py prepare` creates manifests, protocol, and code snapshots.
Submit `run.sh sweep` as array0–20%3; `run.sh select` after the sweep;
`run.sh final` as array0–35%3 after successful selection; `run.sh summary`
after the final array. Selection-seed BANKSY finals reuse selected sweep results.
Failed/missing results block parameter selection or final plotting. Cluster-count
mismatches are retained for inspection, never silently counted as successful.

Outputs: `sweep/`, `selected_parameters.json`, `final/`, `summary/`, `logs/`.
Each fitted run saves labels, embeddings/graph, metrics, graph diagnostics,
and provenance. Summary saves per-seed/mean-SD metrics, paired comparisons to
the earlier batch-zero Harmony results, and a separate PNG/PDF/SVG barplot.
BANKSY comparisons include retuning and do not isolate Harmony's causal effect.
