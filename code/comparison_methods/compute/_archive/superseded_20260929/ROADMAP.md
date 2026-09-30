> Superseded for the active compute benchmark by [PRE_AGGREGATION_PLAN.md](PRE_AGGREGATION_PLAN.md). The implementation now measures native tissue creation before aggregation. This document retains the earlier all-modalities proposal for context.

# Revised compute benchmark: 556 genes

Status: agreed gene count; proposed workflow and execution plan. No revised
benchmark jobs have been submitted. Existing launchers still implement the
historical workflow; do not use them to launch this design unchanged.

## Question and fixed settings

Measure the cost of producing one native 3D simulated tissue and its cell,
16 µm bin, and spot outputs as the number of input cells increases.

- ALBIS, scCube, and SPIDER: exactly 556 generated genes throughout. This is
  already fixed in the historical compute runner, but must also be enforced
  in the revised reference generation, model training, and exported matrices.
  Do not generate 2,000 genes and subset afterward.
- Competitors: fixed 10,000-cell, 556-gene Splatter reference per seed, shared
  between scCube and SPIDER. Record its checksum and include its generation
  cost once in each method's independent end-to-end cost.
- Use current Figure 5A native 3D generation paths and strong-domain-mix
  configuration where applicable, with batch effects disabled. No Harmony,
  clustering, or plotting in the timed simulation pipeline.
- Use the current fixed SPIDER release; pin all package versions and source
  revisions before running. Preserve method-specific expression mechanisms.
- ALBIS/SPIDER: 16 µm bins, circular spots of diameter 55 µm and spacing 100 µm
  where supported by the native path. scCube: native occupancy aggregation
  (3 cells/bin, 10 cells/spot). Record actual capture behavior; these are not
  identical physical measurement models.
- Generate one tissue per method/size/seed, then derive all three modalities.
  Do not regenerate or retrain separately for cell, bin, and spot outputs.
- Keep physical cell density fixed as tissue size grows; record dimensions,
  section thickness/count, and output observation counts. Specify and verify
  the sectioning rule before the pilot so section thickness is not an
  accidental additional scaling variable.
- CPU-only, one computational thread, consistent CPU allocation and hardware
  class where available. Log actual node/CPU and resource limits. Pilot memory
  use determines requests; 64 GB is not a universal capacity ceiling.

## Run order

1. Repair historical accounting: collect array-aware JobID values, preserve
   original adapter errors, and make reconciliation idempotent. Retain old
   outputs separately and label them as the historical workflow.
2. Implement an instrumented runner around the current native 3D paths.
   Separate reference generation, scCube training, tissue generation, and
   slicing/aggregation/export where APIs allow. If an API combines stages,
   report the combined stage rather than inventing separate timings.
3. Pilot 10k and 100k cells, seed 2025, all methods: six tissue runs. Verify
   556 genes at input/training/output, requested cell counts, native XYZ,
   zero batch settings, valid modalities, and complete timing/provenance.
4. Run 10k, 50k, 100k, 200k, 500k, 600k, and 1M cells with seeds 2025, 101,
   and 202. This is 63 tissue runs, each producing all modalities. Valid
   pilot runs can count toward the total if the protocol remains unchanged.
   The 600k point anchors the figure dataset; avoid redundant 20k/40k points.
5. Validate summaries and render runtime/memory scaling plus a stage breakdown
   at 600k cells. Keep failures and resource limits visible in the table.

## Timing and reporting

Primary: from-scratch end-to-end wall time and peak resident memory for one
tissue with all modalities, including reference generation and scCube training.
Train from scratch per independent size/seed benchmark; do not hide a cached
model in the primary runtime. A secondary reuse view may report synthesis cost
with a prepared reference/model, explicitly labeled as such.

Record stage wall time, CPU time, peak RSS, output size, seed, versions,
reference checksum, and scheduler outcome. Sum sequential stage times; take
the maximum stage peak RSS only for sequential stages with no surviving
workers. Record peak memory for the complete process tree if subprocesses
overlap. Include output writing, exclude environment installation and plotting.
For modality-specific totals, include shared tissue cost once; do not add
three modality totals together and triple-count generation.

Show replicate points and mean ± SD for successful runs, with success counts
and failed runs explicit. Treat OOM/timeouts as outcomes under stated budgets,
not observed peak memory/runtime or intrinsic method limits. Three seeds do
not remove hardware variability or establish asymptotic complexity.

## Deliverables

- Versioned protocol/settings and a fresh results directory.
- Pilot validation report and resource recommendations.
- Per-stage logs, timing records, output manifests, and scheduler accounting.
- Tidy measurements CSV; runtime, memory, and stage-breakdown PNG/PDF/SVG.
- Methods paragraph documenting the 556-gene input, reference/training policy,
  native capture differences, replication, and resource limits.
