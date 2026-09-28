# Direct-native compute comparison

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
