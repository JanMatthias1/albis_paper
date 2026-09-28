# Native tissue compute benchmark

The active benchmark measures native tissue creation before aggregation:
ALBIS mRNA instances versus scCube/SPIDER cell-level expression and continuous
3D coordinates. The outputs differ in resolution. See
[PRE_AGGREGATION_PLAN.md](PRE_AGGREGATION_PLAN.md) for the design and feasibility
check. The earlier all-modalities ROADMAP is superseded for this benchmark.

## Protocol

556 generated genes, eight types, fixed cell density, CPU only, one computational
thread, Sapphire Rapids nodes. scCube uses 200 training epochs. Each competitor
uses a shared synthetic 10,000-cell x 556-gene Splatter reference per seed.
ALBIS uses the actual strong-domain-mix expression settings. Full molecule
streams are retained without containment assignment, slicing, or capture.
Native cell expression is materialized before stopping. Output writing and
validation are outside generation timings; compact summaries are saved rather
than full tissue datasets. Molecular count and expression sparsity are recorded.

Primary generation time includes scCube training, labeled as a combined stage.
Input loading/preprocessing, initialization and reference creation are separately
reported; from-scratch totals include these. Memory uses process-tree RSS sampled
at 50 ms plus stage boundaries. Peaks may miss brief transients, and summed RSS
can double-count shared pages. GNU time is retained as a secondary measurement.
Never interpret resource limits as measured maxima or infer complexity from
three seeds. No cached model is used in this benchmark.

## Run sequence

Use a fresh absolute directory; existing runs are never overwritten.

```bash
bash run_compute.sh prepare --root /absolute/new-run
bash run_compute.sh submit --root /absolute/new-run --phase pilot
bash plot_compute.sh --root /absolute/new-run
```

The pilot is six tissue runs: 10k and 100k cells, seed 2025, all three methods.
A reference job precedes competitors. Jobs use two allocated CPUs (one numerical
thread plus monitoring), initially 64 GiB / 24 hours per tissue worker. At most
three tissue jobs run concurrently within a submission, one chain per method.

A dependent review job runs automatically after the pilot. Reporting writes `pilot_review.json` with
conservative resource recommendations. Hardware mismatches or excessive
projected requests hold continuation for review. When all checks pass, the review job submits the remaining grid and schedules final reporting automatically. To resume submission manually after a passed gate:

```bash
bash run_compute.sh submit --root /absolute/new-run --phase full
```

The complete grid is 10k, 50k, 100k, 200k, 500k, 600k, 1M cells, seeds 2025/101/202
(63 tissue runs). Already submitted keys are skipped; reference stages are reused
per seed. A failed chain cancels jobs with invalid dependencies rather than spending
resources on an unchecked larger case; final reporting still runs. Retries require a fresh run directory or explicit
recovery, not silent overwriting. A separate smoke protocol uses 32 output cells,
128 reference cells, one VAE epoch and a reduced SPIDER iteration budget; it must
never be included in the scientific measurements.

## Files and outputs

- `native_worker.py`: native APIs, boundary guards, stage timings and validation.
- `native_benchmark.py`, `native_job.sh`: preparation, reference generation,
  source checks, execution, cluster submission and persistent job records.
- `summarize_native.py`: scheduler reconciliation, CSVs, pilot review and figures.
- `advance_native.py`, `finish_native.sh`: gated continuation and final reporting.
- `run_compute.*`, `plot_compute.*`: entry points for these active scripts.

Run directories include the frozen protocol, implementation audit copies,
settings, per-stage/process measurements, stdout/stderr, source hashes, versions,
reference checksums, scheduler records, `measurements.csv`, `stages.csv`,
`RESULTS.md`, and PNG/PDF/SVG runtime/memory and stage-breakdown plots.
The canonical implementation must match the prepared hashes before a job runs.

## Historical archive

Old code snapshots are in `_archive/historical_20260928`. Old results are in
`sim_paper/data/comparison_methods/compute/_archive/historical_20260928`, including
the legacy `comparison_methods/figure_scaling` tree. Original result paths are
compatibility symlinks so saved provenance continues to resolve. The archive
manifest records every moved tree. These outputs are never read by the new
collector. Old stress helpers are historical and are not current launchers.
