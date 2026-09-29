# Direct-native compute comparison

## Staged 1M, 2M and 5M single-seed run

`native_large_1m_5m_20260929` extends the same protocol using seed 2025 only.
`prepare_large.py` writes scaled inputs and freezes the native workers.
ALBIS/SPIDER are submitted one size at a time after inspecting previous results.
scCube trains once, keeps the VAE in memory, and checks its prior measured peak
before advancing. No model checkpoints are saved. A conservative projection
(prior peak × size ratio × 1.75) must remain below 80% of the allocation.

`sample_rss.py` observes only the scCube worker's `/proc` RSS every 100 ms,
using phase labels written outside timed native calls. This adds a generation
RSS measurement excluding training/validation intervals but including the
resident VAE/reference and allocations retained from earlier stages. It does
not measure incremental allocations or process-tree memory and may miss short
peaks. Cumulative `ru_maxrss`, raw samples, timing stages and gate records are
retained. Read-only finite checks run in chunks to limit validation temporaries.

Regenerate the larger-run table and default three-panel figure with:

```bash
comparison_methods/env/analysis/bin/python sim_paper/code/comparison_methods/compute/report_large.py --root sim_paper/data/comparison_methods/compute/native_large_1m_5m_20260929
```

The plot chooses its range from successful measurements; above 600k it uses a
logarithmic cell-count axis. It distinguishes sampled scCube generation RSS
from earlier cumulative peaks. Missing measurements retain explicit statuses.

## Default figure and retained measurements

The default is the approved three-panel layout:
A, simulation time excluding VAE training and Splatter; B, elapsed time including
Splatter/input preparation and scCube VAE training; C, CPU process peak RSS.
scCube's cumulative memory and separate training batches remain explicitly labeled.

From the project root, regenerate the report and default PNG/PDF/SVG figure with:

```bash
comparison_methods/env/analysis/bin/python sim_paper/code/comparison_methods/compute/report_native_extension.py --root sim_paper/data/comparison_methods/compute/native_extension_200k_600k_20260929
```

The figure is `figures/three_panel_with_setup/compute_three_panel.{png,pdf,svg}`
within that run. Add `--extra-plots` to also regenerate the earlier detailed
runtime, setup, stage, memory and molecule-count figures. Existing extra figures
are retained when the flag is omitted. For plot-only changes, run
`plot_compute_three_panel.py --root <run-directory>` with the same Python.

The full `measurements.csv` and `sccube_setup.csv` remain the data sources for
additional plots: expression generation, spatial placement, native preparation,
input/preprocessing, Splatter, VAE training, first-use totals, process memory,
ALBIS molecule counts, batch identity and source paths are retained. The original
measurement JSON files, validation records, logs and historical script snapshots
remain unchanged. Plotting reads these measurements; it does not run simulations.
Training-time annotations are derived from the setup table rather than fixed text.

## Single-seed extension to 600k (2026-09-29)

`native_extension_200k_600k_20260929` adds 200k, 400k, and 600k cells for
ALBIS, SPIDER and scCube, using seed 2025 only. `prepare_extension.py` scales
the original protocol to these sizes and freezes the executed scripts.
ALBIS/SPIDER use independent native processes. scCube uses the validated split
calls, one new training pass and an in-memory VAE across all three new sizes.
The previous VAE was not saved, so this is a separate training batch.
Large-run output hashing is disabled after the exact 10k equivalence check;
shape and finite-value checks remain. `report_native_extension.py` combines
the earlier pilot with these measurements, retaining batch-specific setup
costs and cumulative scCube memory semantics. Missing measurements are reported.

## scCube split-timing pilot (2026-09-29)

The user approved direct calls to scCube's existing private preparation helper,
`scCube.utils.train_vae`, and `scCube.utils.generate_vae` to separate VAE training
from generation. `sccube_split.py` preserves native outputs and performs no
wrapping, patching, or reimplementation. The original combined-API pilot remains
unchanged. The new run is `sccube_split_pilot_20260929` under the compute data directory.

Two separate processes run the combined API at 10k cells and the split sequence
at 10k/100k. Both use the same 10k-cell reference, 556 genes, seed 2025, 200
epochs, CPU class and one computational thread. The split process trains once,
retains the VAE in memory, and prepares fresh native inputs per size. It saves
no model files. The first generation follows training's RNG state; subsequent
sizes reset generation RNGs to seed 2025. Output digests at 10k must match exactly
before `report_sccube_split.py` accepts measurements and produces plots.

Plots distinguish one-time setup, simulation scaling, and first-use accounting
(setup charged once per size, not independently retrained). The reused Splatter
reference cost comes from its original measurement. scCube RSS is a cumulative
high-water mark over training, earlier generation and validation, not independent
per-size generation memory. Larger sizes and additional seeds remain future work.

## Original combined-API pilot

Each method has its own short script containing top-level calls to its native
public API. There are no custom method functions, adapters, monkey-patching,
replacement implementations, or transformations of returned model outputs.
Input preparation and elapsed-time measurement are permitted by the user.

- `albis_native.py`: `simulate_3d_molecule_sphere_base`. Native output selection
  retains the molecular stream without sectioning or containment aggregation.
- `sccube_native.py`: `pre_process`, `train_vae_and_generate_cell`, and
  `generate_pattern_random`. Training and synthesis are timed together.
  Returned expression and metadata remain unchanged, including native coordinates.
- `spider_native.py`: `simulate_10X_3d` and `get_sim_cell_level_expr`. The returned
  labels, coordinates and expression object, including its native view semantics,
  remain unchanged. Requested cell-type counts are specified as inputs.
- `reference_native.R`: direct Splatter `newSplatParams` and
  `splatSimulateGroups`; native reference counts and Group labels are serialized
  for ordinary input loading. No earlier workflow adapter is called.

The endpoint is native API return, not a standardized output format. No bin/spot
aggregation is performed. Unsupported native outputs remain unavailable: do not
reconstruct fractions, labels, coordinates or other capabilities from other
outputs. No full simulated dataset is exported in this timing experiment; only
read-only shape checks and performance measurements are recorded separately.

## Protocol and interpretation

556 genes; eight input cell types; fixed 10k-cell reference per seed; CPU only;
one computational thread on Sapphire Rapids. scCube uses 200 epochs. ALBIS uses
the strong-domain-mix expression settings. ALBIS/SPIDER tissue dimensions scale
with N. scCube retains its native grid size 8 and coordinate units; we do not
claim matched physical density or rescale its output.

Times cover direct native generation calls, including scCube training. Reference
creation and input preprocessing are reported separately. SPIDER's returned
expression view is not forcibly copied to make its endpoint resemble another
method. This measures the cost of native representations, not equivalent outputs.
Peak RSS is whole-process (including imports/input loading), from the OS; it is
not isolated stage memory or a sum over subprocesses. Scheduler logs can supply
whole-job wall/CPU/memory accounting. Failures are never plotted as successes.

## Execution

```bash
python prepare_native.py --root /absolute/new-directory --smoke
bash submit_native.sh /absolute/new-directory smoke
python report_native.py --root /absolute/new-directory
```

After smoke validation, prepare a fresh directory without `--smoke` and submit
`pilot`. This launches six tissue runs (10k/100k cells, seed 2025, three methods),
plus one reference. Each `.sbatch` executes exactly one method-specific script.
The shell submitter only schedules processes and records job IDs; it never calls
or wraps a model API. Source hashes and snapshots are retained. Do not change
submitted scripts while jobs are pending/running.

The full proposed grid is recorded in `tasks.tsv` (seven sizes, three seeds,
three methods = 63 runs). Remaining points require review of the direct-native
pilot before submission. Smoke settings use 32 output cells, a 128-cell reference,
one scCube epoch and reduced SPIDER iterations; they are not scientific results.

## Excluded implementations

The earlier custom-wrapper code/results are in `_archive/wrapper_based_20260928`
under the code/data compute directories. Their pilot and automatic continuation
were cancelled. Do not reuse these measurements as comparison evidence.
Earlier historical compute results are separately archived.
