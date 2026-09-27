# Figure3 spatial-domain recovery without Harmony

Separate experiment; existing simulations, metrics and figures remain unchanged.
Inputs are the18 existing Figure3 batch-zero QC datasets (strong/weak mix,
cell/bin16um/spot, seeds2025/101/202), read from the existing final_tasks.tsv.

1. prepare.py validates inputs and creates protocol.json, input_audit.json,
   sweep_tasks.json (21 tasks), final_tasks.json (36 tasks), runs.tsv and snapshots.
2. run.sh sweep: BANKSY lambda0,0.1,0.2,0.3,0.5,0.8,1 on strong seed2025.
3. run.sh select: select highest domain ARI among runs with exactly6 clusters;
   exact ties choose smaller lambda. Failure/incompleteness blocks downstream jobs.
4. run.sh final: selected BANKSY and expression-only on both mixes and three seeds.
   Strong BANKSY seed2025 reuses the exact selected sweep result by symlink.
5. run.sh summary: final metrics and revised condensed PNG/PDF/SVG.

BANKSY matches the Figure3 implementation: normalized-to10000, no log, all556
features, max_m1, scaled_gaussian weights, k_geom60/100/8, stagger_scale5 to
separate slices, seeded BANKSY PCA30. Spatial edges are audited for cross-slice
connections. Expression controls use the shared normalize/log/scale/arpack-PCA
helper. Both use 15-neighbour UMAP-weighted expression graphs (NOT Figure5's
K256 graph), fixed random_state0 and the existing Leiden resolution search
(.01..3,15 trials,igraph,n_iterations2) targeting6 domains. No Harmony is called.

Every result saves input/task provenance, embedding-only h5ad, labels.csv.gz,
embeddings_labels.h5ad, graph_diagnostic.json, metrics.json and status.json.
Mismatch runs retain actual ARIs/counts; selection excludes them and the final
figure is withheld if any final run misses6 clusters. No automatic graph fixes.
Parameter selection uses ground truth and seed2025 is not held out.

The revised figure preserves the approved bar style. Left: domain recovery,
new BANKSY/PCA/Leiden and expression/PCA/Leiden, batch0, strong/weak. Right:
existing tuned-batch expression/PCA/Harmony/Leiden cell-type results unchanged.
Mean±sample SD with all three seeds. BANKSY lambda may change during retuning;
previous-Harmony comparisons therefore do not isolate Harmony's effect alone.

Output root: sim_paper/data/figure_3/no_harmony_domain/.
Figure: data/figure_3/ari_recovery_summary/
  banksy_vs_pca_recovery_condensed_no_harmony.{png,pdf,svg}.
From project root, prepare using sim_paper/env/albis-tutorial/bin/python.
Submit run.sh arrays0-20 for sweep and0-35 for final, with selection between;
job IDs/dependencies recorded in jobs.json. run.sh also accepts select/summary.
