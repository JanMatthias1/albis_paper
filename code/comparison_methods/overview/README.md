# Illustrative native 3D outputs at three resolutions

Two scripts: `generate.py` makes the data and calls `plot.py`; `plot.py` can
also redraw saved data. `settings.json` contains the illustrative settings.
This is one fixed realization per method, not a calibrated accuracy benchmark.

Rows: ALBIS, Splatter → scCube, Splatter → Spider. Columns: Cell, Bin, Spot.
Eight cell types, 10,000 cells, 1,200 µm tissue extent and ten sections; each row shares one native 3D
tissue across all three columns. Equal requested type proportions are used for
the competitors. ALBIS uses eight native core/wedge domains with enriched types;
its geometric domain volumes determine the final proportions. Pattern strengths
are illustrative settings, not equivalent numeric targets across algorithms.

## Run

From this directory (or use absolute script paths):

```bash
# One-time dependency download for native circular Spider spots:
bash setup_spider.sh

python generate.py --out /path/to/new/overview_run

# Redraw, without resimulation:
/dcs04/hicks/data/Jan/sim_project/comparison_methods/env/analysis/bin/python plot.py --input /path/to/new/overview_run
```

The parent uses the existing environments in `comparison_methods/env/` and
the existing `splatter_expression.R`. Output goes only to the requested new
directory. It contains method-specific `bin.h5ad`, `spot.h5ad`, `cell.h5ad`,
manifests, settings, the synthetic Splatter pool, logs, and PNG/PDF/SVG figures.
The script refuses to reuse a nonempty output directory. `--method` allows
running an individual method against an already prepared settings/Splatter
directory; that method's output directory must not already exist.

An explicit `--spider-spots square` uses the old installed Spider package and
its working native square aggregation instead; the plot discloses square spots.
There is no automatic fallback and no modification of a native algorithm.

## Native generation and differences to disclose

| Method | Tissue and expression | Bins | Spots |
|---|---|---|---|
| ALBIS | One call to the multiresolution simulator; native default 80 markers/type, 25% shared and 10% noise: **556 genes** | Native molecule-level 16 µm square capture | Native molecule-level circular capture, radius 27.5 µm, spacing 100 µm |
| scCube | VAE trained afresh on 2,000-gene synthetic Splatter input; native organized 3D `generate_pattern_random` | Native `ST` aggregation, target mean 3 cells/unit | Native `Visium` aggregation, target mean 10 cells/unit |
| Spider | Native `simulate_10X_3d`; expression sampled from the synthetic Splatter pool by type | Native 16 µm square aggregation | Native circular aggregation from pinned upstream source; optional installed square mode |

The competitors share the same synthetic input pool, but generate their own
spatial tissue. There are no empirical datasets or empirically pretrained
models. ALBIS's original expression model is retained, not forced to 2,000 genes.

scCube's native lattice has extent eight; its coordinates are multiplied by
the requested tissue extent divided by eight for a declared micrometer scale. Its bin widths are occupancy-driven, not
16 µm bins. Its `Visium` option groups on a square grid and shifts coordinates
to a staggered layout, rather than selecting circular capture regions. Native
scCube expression rescaling is preserved (spot totals >25,000 become 20,000).
These native platform differences are not overridden to make methods match.

For the competitors the wrapper partitions the known native 3D extent into ten
equal z intervals, including the top boundary in the last interval. It calls
native 2D aggregation on each subset independently. Cell z stays continuous;
bin/spot z is the physical section midpoint. `slice_id` is zero-based. No
independent 2D layouts are aligned or invented, and no misalignment is applied.

Colors show cell type; mixed capture units show dominant type. ALBIS dominance
uses native molecule fractions; comparator dominance uses contributing cell
counts. Empty units are saved but omitted from the figure. All panels use one
viewing angle and coordinate frame; native scCube coordinates and ALBIS's
original centered coordinates are also preserved. Marker sizes are illustrative,
not physical capture footprints.

## Spider version

The installed 1.1.0 circular branch failed a direct four-cell test in its
neighbor helper (`n_neighs=None`); source inspection also found the incorrect
membership slice `sn[Num_sample:, Num_sample]`. Upstream's **sim_expr module**
corrects this to `sn[Num_sample:, :Num_sample]` and supplies a revised neighbor
helper. The upstream top-level alias still points to the old sim_naive branch,
so the wrapper explicitly imports `spider.sim_expr.get_sim_spot_level_expr`.
`setup_spider.sh` obtains the **unmodified complete Python package sources** at
commit `6ccd4da77257f2807c430f8f42fbe2dc175991de` in
`sim_paper/env/spider-overview-src`. This isolates the newer implementation from
the old benchmark. Manifests identify the imported path and upstream commit.

Source: https://github.com/YANG-ERA/Spider/tree/6ccd4da77257f2807c430f8f42fbe2dc175991de

## Checks

Generation checks finite nonnegative expression, 3D coordinates, all requested
sections, expected labels and unique observation IDs. Spider square aggregation
must assign each cell exactly once and conserve gene totals; scCube membership
must cover the source section without duplicated cells. All three outputs are
derived in one method run; manifests report the realized gene/observation counts.
Circular Spider membership and expression are checked independently against
Euclidean distance to spot centers and the sum over cells inside each circle.
Pilot settings can be supplied using `--settings`; do not interpret a one-epoch
VAE or short Spider optimization smoke test as a manuscript-quality realization.

Verified 2026-09-21: an end-to-end run with **800 cells, eight types, ten sections,
one VAE epoch and 100 Spider iterations** completed all nine outputs, native
circular membership/expression checks, and plotting. ALBIS retained 556 genes;
the smoke-test Splatter pool used 80 genes. Outputs and all logs are under
`sim_paper/data/comparison_methods/overview_smoke_20260921/`. The default
10,000-cell/2,000-Splatter-gene run completed as job 35831924 at 600 µm extent. The smoke preview is
labeled using `plot.py --label 'Smoke test ...'`.

The tissue extent was subsequently doubled to 1,200 µm at the user's request.
Cell count, ten-section count, bin width and spot capture sizes are unchanged.
This increases physical spot positions for ALBIS/Spider but lowers cell density;
scCube's occupancy-driven native spot count remains approximately unchanged.
All sections are plotted together. Saved cells retain continuous z; saved capture
units use physical section midpoints. The updated **display** projects cells onto
their section planes at their physical section spacing, so all three resolutions
visibly show the same ten-section stack.
Stored section IDs remain 0–9. The display follows `applications_albis/alignment/plot_figure4a.py`:
white background, no axes or grid, small translucent points and one shared cell-type legend.
This is a visualization transformation only, disclosed beneath the figure.
Use `--section-spacing` to change the display gap and `--output-prefix` to save
an alternative rendering without overwriting the existing overview files.
Existing datasets retain their own settings and are not rewritten.

The original 600 µm full run was redrawn as `overview_sections.{png,pdf,svg}`
with the new Cell → Bin → Spot order. The 1,200 µm rerun has not been submitted
successfully; the redraw therefore still shows the existing 600 µm dataset.
The compact 12 × 4.3 in rendering, matching the Figure 4a panel dimensions, is
saved as `overview_sections_compact_v10.{png,pdf,svg}`. Its display projection
compresses the z-axis modestly to fit the three section stacks in the panel;
this does not alter the saved coordinates.
Version 10 removes the square Axes3D clipping patch from the scatter artists
so the full tissue silhouettes remain visible, and uses two rows of four cell
types in the legend. To reproduce it from the repository root:

```bash
comparison_methods/env/analysis/bin/python sim_paper/code/comparison_methods/overview/plot.py \
  --input sim_paper/data/comparison_methods/overview_full_35831924 \
  --output-prefix overview_sections_compact_v10
```

Version 11 increases the display section spacing by 25% to make the slices
easier to distinguish, retaining the unclipped silhouettes and two-row legend:

```bash
comparison_methods/env/analysis/bin/python sim_paper/code/comparison_methods/overview/plot.py \
  --input sim_paper/data/comparison_methods/overview_full_35831924 \
  --section-spacing 1.25 --output-prefix overview_sections_compact_v11
```

This changes the display z scale only; the stored coordinates are unchanged.

Version 12 also trims the unused left and right canvas margins, retaining a
small white border and the 25% larger display section spacing:

```bash
comparison_methods/env/analysis/bin/python sim_paper/code/comparison_methods/overview/plot.py \
  --input sim_paper/data/comparison_methods/overview_full_35831924 \
  --section-spacing 1.25 --trim-side-margins --output-prefix overview_sections_compact_v12
```

To show the fifth section as a top-down 2D view in the same nine-panel layout:

```bash
comparison_methods/env/analysis/bin/python sim_paper/code/comparison_methods/overview/plot.py \
  --input sim_paper/data/comparison_methods/overview_full_35831924 \
  --slice-number 5 --trim-side-margins --output-prefix overview_slice5_compact
```

`--slice-number` is one-based: slice 5 selects stored `slice_id=4` for every
method and modality. This redraw uses the existing data and shared x/y limits.
