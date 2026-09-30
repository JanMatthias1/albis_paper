# Prior testing: ALBIS vs. Spider vs. scCube (from `/dcs04/hicks/data/Jan/sim_project/comparison_methods`)

Working notes, not manuscript content. Written 2026-09-21 to summarize what was
already tested in the old scratch repo (`comparison_methods/`, no git, separate
from `sim_paper/`) before starting fresh work here on the new Figure 4A/4B
(ALBIS vs. scCube table + cross-technology slice comparison). Everything below
is either read directly from that repo's code/data or from
`project_figure5_comparison_methods` session memory — verified against the
actual adapter source where it mattered, not assumed from either method's
public docs (the adapters' own docstrings say the same: the public scCube/Spider
READMEs implied APIs that don't actually exist in the installed packages).

Full design doc, if you want the complete methodology: `comparison_methods/README_FIGURE5.md`.

## What this was: Figure 5, not Figure 4

The old repo built an end-to-end **reference-free** benchmark ("Figure 5")
comparing three complete workflows — `our_method` (ALBIS), `Splatter→Spider`,
and `Splatter→scCube` — on how well each hits user-*requested* spatial/molecular
targets (no empirical ground truth involved; everything scored against the
requested specification, not biological realism). It reached a finished draft:
Stage 1 calibration (20 candidates × 5 seeds per scenario/org-level, all 3
workflows) and Stage 2 held-out evaluation (20 replicates, technology=cell
only) both completed 2026-09-07, with figures 5a/5b/5c drafted.

Your new Figure 4A/4B is a **different, narrower ask** — just ALBIS vs. scCube,
starting with a comparison table (4A) and a slice-level qualitative comparison
across technologies (4B), not the full controllability/adherence benchmark.
But the same two methods are involved, so every hard technical limitation
found below still applies directly.

## The one limitation that will shape 4B most: scCube is 2D-only outside `random_null`

This is the single most load-bearing finding for what you're about to build,
and it's verified against scCube's actual installed source
(`env/scCube_src/scCube/sccube.py`), not assumed:

- `scCube.generate_pattern_random(..., spatial_dim=3)` is genuinely 3D — but
  it's only used for the **random_null** scenario (no intended cell-type
  organization).
- The **shaped** generators — `generate_pattern_custom_cluster` / `_ring` /
  `_stripes` (used for clustered / layered / mixed_interface, i.e. any
  scenario with actual spatial structure) — **hardcode 2D**. Confirmed by
  reading `__generate_custom_helper`: it only ever writes `point_x`/`point_y`,
  there is no z-column path at all, not even an unused one.
- `generate_spot_data` (scCube's own bin/spot aggregation) **hard-asserts
  `spatial_dim == 2`** — so scCube cannot produce 3D bin/spot data under any
  scenario, including random_null.

Net effect: **scCube is only ever 3D for one specific case (random_null,
cell-level).** Every other technology×scenario combination it produces is 2D.
The old adapter records this honestly per-replicate
(`spatial_dimensions_supported` in the manifest) rather than faking a z-axis —
worth doing the same if you extend this.

**Spider has the identical restriction at bin/spot** (verified against
`env/spider/lib/.../spider/` source): `spider.simulate_10X_3d` is real 3D for
`cell` technology, but `spider.get_sim_spot_level_expr`'s square-grid
aggregation is 2D-only — "same restriction scCube's `generate_spot_data` has"
per the adapter's own comment.

**ALBIS (`our_method`) was the only one of the three that was fully 3D at
every resolution** (cell, bin, spot) in the old benchmark.

**Implication for your 4B ("plot of the slice... across the different
technologies"):** if the panel is a genuine 2D *slice* (not a full 3D render),
scCube can actually be compared at every technology and scenario — a slice
sidesteps the 3D-vs-2D asymmetry. If instead you want a 3D structural
comparison, scCube will only be an honest comparison for the random_null
scenario at cell resolution; anywhere else you'd be comparing ALBIS's real 3D
structure against scCube's 2D output, which needs to be labeled as such (the
old Figure 5a panel literally tags those panels "2D output" rather than hiding
the asymmetry — see below).

## scCube integration was genuinely finicky — two real adapter bugs found and fixed

Not scCube being slow to use correctly — actual wrong-output bugs, caught only
because Stage 1 calibration numbers looked off and someone dug into it (see
`code/slurm/rerun_sccube_shaped_array.sh` header):

- `generate_pattern_custom_stripes` (layered scenario) was called without a
  required `infiltration_*_list` argument in the first integration pass —
  silently produced wrong structure until fixed.
- `generate_pattern_custom_ring` (mixed_interface scenario) was called with a
  flat `ring_celltype_list` and a dead infiltration knob that scCube's real
  API doesn't use that way.

Both fixed in `code/workflows/sccube_adapter.py`, then a targeted rerun
(`rerun_sccube_shaped_array.sh`, `aggregate_sccube_rerun.sh`) redid just the 6
affected calibration combos. **Lesson for reuse:** don't trust a first-pass
scCube shaped-pattern integration without checking the achieved spatial
statistics against the request — the API accepts wrong arguments silently
rather than erroring.

## scCube runtime is slow, highly variable, and prone to cluster hangs

From `figure5/metrics/summary/computational_burden.csv` (technology=cell,
20 replicates per combo):

| workflow | median runtime | range seen | peak memory |
|---|---|---|---|
| our_method (ALBIS) | ~90–250 s | tight | ~5.9 GB |
| Splatter→Spider | ~90–400 s | tight | ~1.2 GB |
| Splatter→scCube | **1,200–11,600 s** (20 min – 3.2 h) | **CI tail up to 32,530 s (~9 h)** on random_null | ~1.3 GB |

scCube trains a VAE from scratch per replicate (confirmed: no pretrained
weights are touched in the reference-free path — `train_vae_and_generate_cell`
trains fresh every time; the ~300-model pretrained atlas is a separate, unused
code path). That VAE training is both slow on average and has a long,
unpredictable tail.

**Separately, a real cluster gotcha, not a scCube algorithmic issue:** 7 of 70
Stage-2 scCube evaluation tasks hung for 17+ hours (`TotalCPU 00:00:00`, node/NFS
stall) from multiple VAE jobs landing on the same bad compute nodes
(`compute-111`, `compute-095`, `compute-116`). `--cpus-per-task` did not
prevent SLURM from co-locating them; only `--exclude=<those nodes>` fixed it.
**If you submit any new scCube array here, exclude those three node names
preemptively** rather than rediscovering the hang.

## Controllability: ALBIS tracks requests closely; scCube is inconsistent; Spider under-responds

From Stage 2 held-out evaluation (calibration slope = 1.0 means perfect
tracking of requested vs. achieved organization strength):

- **our_method (ALBIS):** slope ~1.0–1.03, r ~0.99 on the primary statistic
  (same-type neighbor enrichment) for clustered and layered. Tracks requests
  closely across the full low/medium/high range.
- **Splatter→scCube:** **inconsistent by scenario** — slope 0.27 for
  clustered, but slope 1.00 for layered. No simple rule for when it will vs.
  won't track a request; each scenario needs its own check, don't assume
  behavior transfers.
- **Splatter→Spider:** systematically **under-responds** — slope ~0.40,
  plateaus around 0.49 regardless of how high the request goes. Also a
  Stage-1 finding: both Spider and scCube hit a **controllability ceiling at
  "high" organization** (e.g., spider mixed_interface_high combined_error
  0.52) — they fit low/medium fine but can't reach the top of the requested
  range. ALBIS did not show this ceiling.

## Universal, unresolved limitation: sequencing depth is unwired everywhere

None of the three workflows honored the requested counts/cell (500 requested):
ALBIS realized ~8,000, Spider ~60,000, scCube ~1,700 — all off by different,
large multiples, in different directions. This was never wired into any of
the three adapters (`sequencing_depth` parameter is a no-op end to end).
**If depth realism matters for your new comparison, this needs fixing in
whichever adapter(s) you reuse — it's not fixed anywhere currently.**

## Scope gaps: what was never actually done

- **Only `technology=cell` got a full independent Stage-2 evaluation.**
  bin/spot were deliberately **not** separately calibrated or evaluated —
  the plan was to reuse the same frozen cell-level spatial configs on the
  (stated, not verified) assumption that "spatial knobs are technology-
  independent." That assumption was never tested end-to-end for bin/spot.
- **Cross-technology representative renders were flagged as needed but never
  built.** The existing `figure5a.py` renders 2 scenarios × 3 *workflows* at
  one fixed technology (cell) — it does **not** compare technologies against
  each other. The old session's own remaining-work note says explicitly:
  "bin/spot representative renders for 5a (need ~6 gens)" — this is almost
  exactly your new 4B ask, and it was never done. There's no existing bin/spot
  generation output to reuse for this; it would need fresh runs.

## Reusable assets, if you want a starting point rather than building from scratch

- `code/workflows/our_method.py`, `sccube_adapter.py`, `spider_adapter.py` —
  working adapters with a shared contract (`code/common/contract.py`) taking
  the same high-level spec (n_cells, n_genes, cell-type proportions, scenario,
  organization level, technology, seed) and producing a common output schema.
  scCube's and Spider's adapters both record `spatial_dimensions_supported`
  per replicate — reuse that field rather than re-deriving the 2D/3D fact
  above.
- `code/figures/figure5a.py` — already builds a multi-panel representative-
  output comparison (per-axis min-max normalized, since the three engines use
  incompatible native coordinate systems: ALBIS microns centered at 0, Spider
  microns from a corner, scCube an 8³ lattice) with 2D-output panels
  explicitly tagged. Good structural template for a slice comparison, but
  it compares *workflows* at fixed technology, not *technologies* — would
  need adapting, not just rerunning.
- `figure5/evaluation/<workflow>/<scenario>_<org>_cell/rep0/output.h5ad` —
  existing generated outputs at cell resolution, if a quick illustrative
  slice at that resolution is enough to start with before generating anything
  new for bin/spot.
- Environments: `env/our_method`, `env/sccube`, `env/spider`, `env/analysis`
  (plotting/metrics), `env/splatter` (R, expression input for Spider/scCube).

## Suggested questions to settle before building 4B

1. **2D slice or 3D render?** Determines whether scCube can be shown honestly
   across all technologies/scenarios (slice) or only for random_null at cell
   resolution (3D render).
2. **Which technologies?** bin/spot have zero prior independent evaluation for
   either method here — this would be new work, not a rerun.
3. **Which scenario(s)?** Given scCube's inconsistency, picking a scenario
   where it happens to track well (layered) vs. poorly (clustered) will tell
   very different stories — worth being deliberate about which one(s) 4B
   shows, and saying why in the caption.
4. **Runtime budget** — if scCube is involved and you're generating fresh data
   (not reusing rep0 outputs above), budget for the multi-hour tail and
   exclude the three known-bad compute nodes from the start.
