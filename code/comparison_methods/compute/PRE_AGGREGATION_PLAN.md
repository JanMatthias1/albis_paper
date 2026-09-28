> Updated execution constraint (2026-09-28): the active implementation is
> described in README.md and uses direct native API calls only. Earlier suggestions
> below about instrumentation wrappers, forced expression materialization, or
> standardized physical coordinates are superseded. Native outputs remain
> unchanged; missing capabilities are unavailable. Wrapper-based results are excluded.

# Native tissue generation before aggregation

Status: implementation and three-method smoke validation completed on 2026-09-28.
The active run is recorded in sim_paper/data/comparison_methods/compute/CURRENT_RUN.txt. This plan narrows the older ROADMAP.md proposal:
cell/bin/spot exports and sectioning are outside the measurement endpoint.

## Scientific question and endpoints

How much time and memory does each method require to construct its native
3D tissue representation as requested tissue cell count increases?

- ALBIS: fully materialized molecule XYZ, gene IDs, source type/domain labels,
  plus the native tissue cell metadata. Include tissue placement, expression
  count generation, molecule placement, and assembly of returned arrays.
- scCube: fully materialized cell expression, cell types, and continuous XYZ.
  Include native expression synthesis and 3D spatial pattern generation.
- SPIDER: fully materialized cell expression, cell types, and continuous XYZ.
  Include native 3D placement/type assignment and reference expression sampling.

These outputs have different resolution. Report this as the computational cost
of native representations, not a speed ranking for identical output. Do not
call cell-expression counts explicit molecules or divide all runtimes by a
common molecule count. Requested cells is the shared horizontal axis.
Exclude slicing, molecular containment assignment, observed-cell aggregation,
bin/spot capture, batch effects, downstream analysis, plotting, and final export.
Include necessary in-memory output construction. Store compact validation and
measurement records outside the timed generation region.

## Feasibility findings

### ALBIS: existing API supports the endpoint

`albis/albis/simulation_sphere.py:simulate_3d_molecule_sphere_base` returns before
sectioning. Its default output selection still performs containment assignment
and aggregation to observed cells, so the default call is inappropriate here.
Calling the base API with `output_modalities=("bin", "spot")` retains the full
molecule stream and skips that cell work. These flags select retained streams;
they do not perform bin/spot aggregation when using the base API.
The base still constructs native true-count/placeholder AnnData objects; include
that existing API overhead and do not label the placeholder observed counts as
an observed-cell simulation. No simulator source modification is required.

Executed a 32-cell smoke check using the strong-domain-mix configuration:
556 genes and 2,828 molecules. All four full-stream arrays matched the default
base call exactly. Replacing containment and all three aggregation functions
with failure sentinels confirmed none was called in the molecular-only path.
This validates the API boundary, not large-run runtime or memory feasibility.

### SPIDER: separate native calls

The installed st-spider 1.2.0 path used in overview/generate.py calls
`simulate_10X_3d`, then `get_sim_cell_level_expr`. Stop after these and explicitly
materialize the returned expression: the latter returns an AnnData view, so
stopping at view creation would omit part of the required work. No capture call
is needed. Its >10,000-cell spatial path differs from the small-data path, so
pilot both sides of this threshold. Retain native optimization settings and
record convergence diagnostics where available.

### scCube: native generation and unsliced placement are available

Use `train_vae_and_generate_cell`, then `generate_pattern_random` with
`spatial_dim=3, is_split=False`. Stop before `cell_data`/modality conversion and
`generate_spot_data_random`. Training and synthesis are combined in the public
API; the installed implementation calls separate `train_vae` and `generate_vae`
functions. Start with a combined training+synthesis measurement and a separate
spatial-placement measurement. Optional transparent timing wrappers can split
those internal calls after an equivalence check; do not rewrite the algorithm.
`load_vae_and_generate_cell` also supports a separately labeled reuse experiment.
Never estimate reuse cost by subtracting training time from another run.

Competitor feasibility is based on installed-source inspection and the existing
native workflow, not new benchmark executions. Dense scCube outputs and ALBIS
molecular arrays make larger-run memory requirements a pilot question.

## Proposed controlled settings

- 556 genes generated from the start; eight cell types; fixed method-specific
  expression and spatial settings across sizes; zero batch perturbations.
- Fresh synthetic Splatter reference: 10,000 cells x 556 genes per seed, shared
  by scCube/SPIDER. Verify dimensions and checksum. The current Figure 5A pool
  is 10,000 x 2,000 and cannot be reused as the matched benchmark input.
- ALBIS: use the actual strong-domain-mix source configuration, including its
  expression depth and cell radii, rather than the overview helper defaults.
  The expression settings control molecule count and therefore workload.
- scCube: retain 200 training epochs and the current native pattern parameters;
  report these as a fixed benchmark configuration, not a convergence guarantee.
- Preserve native sphere/cube geometries. Proposed common physical density:
  rho = 600000 / ((4*pi/3)*2050**3). At N cells set ALBIS radius to
  (3*N/(4*pi*rho))**(1/3) and competitor cube extent to (N/rho)**(1/3).
  This deliberately differs from giving a sphere and cube the same bounding
  width, which would not match density. Keep cell radii fixed; scale domain
  length parameters with tissue extent. Record scCube native-coordinate scaling
  and its fixed latent grid/pattern parameters. Validate these choices in pilot.
- CPU only, one computational thread, consistent CPU class and allocation;
  record actual node, CPU, thread settings, package versions and source hashes.
- Use fresh processes, a documented uniform cache policy, and no model reuse
  in the from-scratch benchmark. Record imports/loading separately from the
  generation interval, including any lazy initialization inside native calls.

## Measurements and accounting

Report native generation wall time (including scCube training, explicitly
labeled), stage CPU time, and process-tree peak resident memory through the
pre-aggregation endpoint. Keep model training visible in the stage breakdown.
Record reference generation separately and also show the from-scratch total
including it once for each competitor. ALBIS requires no external reference.
Reference cost can be shared physically while being charged once to each
independent method total. Exclude queue time and installation.

Measure the whole worker and stage boundaries with monotonic timers and a
process-tree memory monitor. Avoid inferring stage memory by subtracting
cumulative high-water marks. Sum sequential times; use maximum stage peaks only
when stages are separate processes without surviving workers. Report memory
sampling frequency and any measurement limitations.

Record requested/actual cells, genes, molecule count for ALBIS, expression
nonzeros/dtypes, dimensions, seed, reference checksum, successful endpoint,
resource limits, and scheduler outcome. Validate finite XYZ and expression,
nonnegative expression, row alignment, all genes, and no aggregation calls.
Keep failed/OOM/timed-out runs visible; limits are not measured capacity.

## Execution plan

1. Implement a dedicated pre-aggregation runner in this compute directory,
   using the existing isolated environments. Leave historical runners/results
   intact. Persist the resolved protocol and planned stages before execution.
2. Add small boundary/equivalence checks: ALBIS full-stream parity; materialized
   SPIDER expression; scCube expression/coordinate row alignment; guards that
   sectioning/capture functions are never called. Validate any timing wrappers.
3. Pilot 10k and 100k cells, seed 2025, all methods: six tissue runs plus one
   10k x 556 reference. Use allocated compute resources. Verify endpoints,
   runtime accounting, gene counts, spatial settings and memory; choose larger
   resource requests from measured usage. Do not assume a 64 GiB ceiling.
4. Freeze the protocol after pilot review. Proposed main grid: 10k, 50k, 100k,
   200k, 500k, 600k, 1M cells; seeds 2025, 101, 202 = 63 tissue runs and three
   reference preparations. Valid unchanged pilot runs may count. Stage larger
   runs according to observed resource use rather than launching all blindly.
5. Export tidy measurements and runtime/memory scaling plots, with individual
   replicate points and success counts. Add a stage breakdown at 600k and an
   ALBIS molecule-count column/panel. Clearly label molecular versus cell-level
   outputs and include reference/training costs in captions and methods text.

First implementation milestone: the instrumented runner and six-run pilot,
not the full sweep. Large-scale feasibility and a reliable compute budget
remain unmeasured until that pilot finishes.
