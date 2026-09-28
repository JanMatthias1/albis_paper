# Preliminary compute comparison

Run commands from this directory. The Python runner uses the existing workflow
adapters and isolated environments in the project's `comparison_methods/` tree;
these files are not a standalone installation.

## Plot completed results

```bash
bash plot_compute.sh
```

Reads the corrected `comparison_methods/figure_scaling/raw` sweep, never the
archived fixed-volume results. Writes PNG, PDF, SVG and a measurement CSV to
`sim_paper/data/comparison_methods/compute/figures`.

```bash
bash plot_compute.sh --input-dir /path/to/raw --output-dir /path/to/figures
```

Columns are cell/bin/spot; rows show wall time (log scale) and peak resident
memory (GiB), with common y limits across technologies. SPIDER and scCube include
Splatter: sequential wall/CPU times are summed; peak memory is the maximum of
stage peaks, not their sum. Failure crosses sit outside the metric axes; failed
runs are excluded from curves and retained in CSV. CSV also includes CPU seconds,
source paths, seeds and failure reasons. One run per point: no uncertainty bands
or fitted complexity claims. Lines only connect observations. Input size is
requested cells, not output spots/bins. Some competitor spatial paths are 2D;
this is a comparison of the configured workflows, not identical 3D algorithms.

## Run compute

One Python file dispatches all three methods (or selected methods) and runs
Splatter once when needed. The default design matches the completed sweep:
556 genes, random_null, CPU, deterministic seed, density-preserving geometry.
scCube trains its VAE from scratch. Override `--seed` for a new replicate and
use a separate output directory.

```bash
# Preview all subprocess commands without running or writing outputs.
bash run_compute.sh --n-cells 10000 --technology cell --dry-run

# One point, all methods (run on allocated compute resources).
bash run_compute.sh --n-cells 10000 --technology cell

# Select a method and a fresh output directory.
bash run_compute.sh --n-cells 20000 --technology spot --method albis --out-dir /path/to/new-point

# Submit all 15 size/technology points, all methods.
sbatch run_compute.sh

# A selected method across the sweep, in a separate output tree.
COMPUTE_OUTPUT_ROOT=/path/to/new-sweep sbatch run_compute.sh --method albis
```

Default new results: `sim_paper/data/comparison_methods/compute/raw`.
Existing nonempty point directories are refused. Each step retains stdout,
stderr, GNU time measurements and its workflow manifest. The summary is saved
after every step, including failures. A failed Splatter stage skips dependent
methods while ALBIS can still run. The runner exits nonzero if any step fails.
SPIDER's known upstream 3D reshape failure above 10k cells is expected; the
launcher provides the full grid for reproducibility, not an upstream fix.
Slurm logs appear in the submission directory; resources match the prior sweep
(4 CPUs, 64 GiB request, two-day limit). No full sweep was submitted when creating
these scripts.

## Larger stress sweep

`prepare_stress.py`, `run_stress.sh`, and `finish_stress.sh` extend the sweep to
200,000, 500,000 and 1,000,000 cells, for cell/bin/spot and all three methods.
There are 27 independent method jobs, with at most three running concurrently.
Each gets 64 GB, four allocated CPUs and 48 hours; computational thread counts
are fixed to one. The wall limit applies to each method pipeline, including
its Splatter stage, rather than all three methods sharing one job's limit.
Each competitor generates a fresh full-size Splatter pool; its costs remain
included. This differs from the illustration's fixed 10,000-cell reference.

From the project root, choose a fresh absolute output root:

```bash
root=/dcs04/hicks/data/Jan/sim_project/sim_paper/data/comparison_methods/compute/stress_YYYYMMDD
comparison_methods/env/analysis/bin/python sim_paper/code/comparison_methods/compute/prepare_stress.py --root "$root"
sbatch --output="$root/logs/%A_%a.out" sim_paper/code/comparison_methods/compute/run_stress.sh "$root"
# Replace ARRAY_JOB_ID with the returned ID. Runs after success OR failure.
sbatch --dependency=afterany:ARRAY_JOB_ID --output="$root/logs/finalize_%j.out" sim_paper/code/comparison_methods/compute/finish_stress.sh "$root" ARRAY_JOB_ID
```

The new root links the original 15 baseline directories without modifying them.
The dependent job saves Slurm accounting, classifies unfinished method stages
(timeout, out of memory, cancellation, infrastructure failure, or incomplete),
and writes combined baseline + stress plots and CSV under `figures/`.
Adapter-reported errors retain their failure reason. Missing measurements are
never treated as successful runs or substituted with the requested memory limit.
SPIDER cell tests above 10k are expected to fail in its native 3D path; the
bin/spot paths are separate 2D workflows. These are workflow stress tests,
not matched native 3D simulations. One seed per point is exploratory evidence,
not a precise maximum capacity estimate.

## Capture correction: 16 µm bins and 55 µm circles

The current compute runner overrides the historical benchmark contract to use
16 × 16 µm bins and 55 µm diameter circular spots at 100 µm center spacing.
ALBIS consumes these dimensions directly. For SPIDER spot aggregation only,
the adapter loads the unmodified upstream `sim_expr` implementation at commit
`6ccd4da77257f2807c430f8f42fbe2dc175991de` under an isolated module name.
Installed SPIDER placement stays unchanged; contracts without the explicit
`pinned_circle_v1` flag retain legacy square behavior. Capture geometry and
source commit are stored in the output AnnData; the resolved contract is in
its manifest. scCube still uses its native occupancy-driven aggregation
(three cells/bin, ten cells/spot); physical dimensions are not enforced there.

The replacement sweep is now **72 independent tests**: eight sizes
(10k, 20k, 40k, 50k, 100k, 200k, 500k, 1M), three technologies, three methods.
`prepare_stress.py` no longer links the old baseline, because its capture
geometry differs. The old `stress_20260923` jobs 35882445/35882478 were cancelled;
their partial results remain archived. Use a fresh root such as
`stress_bin16_circle55_20260923` with the submission commands above.
Summaries and CSV record a geometry version; the plotter refuses mixed versions.
The earlier stress-sweep section describes the initial submission; these
geometry and baseline-rerun changes supersede its 27-job/linking description.
