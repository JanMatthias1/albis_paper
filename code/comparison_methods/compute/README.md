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
