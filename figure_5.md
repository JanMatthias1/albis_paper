# Figure 5 — comparison between simulation methods

Planning notes saved 2026-09-26 for the next review.

Update 2026-09-27: installed `st-spider==1.2.0` in the SPIDER environment and
removed the old pinned-source override from the overview generator. Native
3D generation at 20,000 cells passed a one-iteration smoke test (finite 3D
coordinates and correct label count). The earlier large-3D failures below
describe the previous implementation; the full 600,000-cell overview subsequently completed with the release fix
in job 35984840; all nine panels are verified.

## Guiding principle: honest and fair comparison

The objective is to characterize each method accurately, including strengths,
limitations and native design choices. A result favoring a comparator is a
valid result, not a reason to alter the benchmark. Do not design the analysis
to favor ALBIS or interpret recoverability as a proxy for realism.

- Separate **native-capability comparisons** from **controlled sensitivity
  analyses**. Native runs retain each method's intended behavior and disclose
  differences; controlled runs match specified factors and disclose changes.
  Neither should be presented as controlling factors it does not control.
- Define the scientific question, inputs, primary metrics and tuning rule
  before evaluating the new comparison. Document that initial pilots have
  already been inspected; subsequent protocol choices are informed by them.
- Use the same downstream implementation and tuning budget where applicable.
  Known-eight-type Leiden tuning uses ground-truth category count and must be
  labeled accordingly. Do not choose resolutions, genes, seeds or panels to
  maximize a preferred method's reported performance.
- Do not force an unsuitable preprocessing step solely for procedural symmetry.
  For these runs without an injected technical batch effect, report paired
  PCA-only and PCA-plus-Harmony analyses for every method on matched inputs.
  Treat batch-correction performance as a separate controlled experiment.
- Match cell count, input gene count, reference protocol and evaluation target
  where feasible. Report remaining differences in expression depth, detection,
  marker strength, capture geometry and occupancy. Equal gene count alone
  does not make the datasets equally difficult. Synthetic feature names are
  not shared biological gene identities across ALBIS and Splatter.
- Resolve cell-fraction versus molecule-fraction scoring before a direct
  deconvolution ranking. Report coverage and filtering, including failures,
  alongside accuracy; high accuracy on fewer retained spots is incomplete.
- Repeat simulations using a predefined set of independent seeds and report
  variation. Distinguish simulation replicates from repeated gene subsets or
  clustering seeds; spots within one simulation are not independent simulation
  replicates. The current single-run pilots remain exploratory.
- Report runtime and memory with hardware, parallelism and stage boundaries.
  Include required training and preparation, or show them separately. Record
  unsuccessful runs and resource limits without claiming they prove a general
  limitation of a method.
- Keep source versions, configurations, seeds, input provenance and complete
  results. Show representative panels using a declared selection rule and
  retain all results in supporting material. Clearly label custom transforms,
  rounding, or other adaptations and preserve native baselines.

This is the intended protocol for further work, not a claim that the existing
pilots already meet all these controls. The final paper should state what can
and cannot be inferred from each panel.

## Proposed panels

1. **Visual comparison between methods** — ALBIS, scCube, and SPIDER at cell,
   bin, and spot resolutions, with the section stack and one section view.
   Canonical 600,000-cell comparison: `data/comparison_methods/figure_5A_600k/`.
   The current ALBIS panels use **uncropped**
   `figure_4/cross_modality_alignment/strong_domain_mix_shift3x/data`, with
   aligned `spatial_3d` coordinates and full-field 16 µm bins. The previous
   cropped rendering is preserved in
   `_archive/figure_5A_cropped_before_uncropped_20260927/` under comparison_methods.
   SPIDER regeneration with st-spider 1.2.0 and rendering completed as job
   35984840 (600,000 cells, 660,490 bins, 16,810 spots before empty-unit filtering). ALBIS and scCube reuse the validated 600k sources below. Prior
   previews and their plotting summaries were moved into
   `data/comparison_methods/_archive/figure_5A_pre_600k_all_methods_20260927/`;
   the archive has a complete move manifest. Source datasets stay in place.
   Latest three-method rendering (2026-09-27):
   [overview with native 10k SPIDER](data/comparison_methods/_archive/figure_5A_pre_600k_all_methods_20260927/overview_shared600k_spider10k_20260927/overview_combined.png).
   This uses the same ALBIS/scCube inputs described below and the original
   `overview_full_35831924/spider` 10,000-cell native 3D tissue for **all three**
   SPIDER modalities. It does not stack independently generated 2D tissues.
   Row labels disclose cell counts. SPIDER uses its native 600 µm extent,
   enlarged within the panels; ALBIS/scCube retain 4,100 µm extents. The footer
   discloses the different scales. Section 5 is stored ID 4 for every row.
   Updated reference (2026-09-27):
   [overview_combined.png](data/comparison_methods/_archive/figure_5A_pre_600k_all_methods_20260927/overview_shared600k_20260927/overview_combined.png).
   ALBIS and scCube each contain **600,000 cells total across ten sections**.
   ALBIS uses the exact shared-tissue raw cell/bin16um/spot files from
   `figure_4/cross_modality_alignment/strong_domain_mix_shift3x_cropped_bin/data`,
   displayed with aligned `spatial_3d` coordinates. Its 16 µm bins retain the
   2,221 µm crop; cells and spots retain their original fields. scCube reuses
   the completed `overview_figure2_600k_35881874/sccube` realization, with
   188,896 bins and 59,285 spots at targets of 3 cells/bin and 10 cells/spot.
   These are occupancy-driven units, not matched physical capture footprints.
   Both rows use the same physical display scale and section 5 (stored ID 4).
   Genes and proportions are not matched: ALBIS has 556 genes and native
   realized proportions; scCube has 2,000 genes and equal requested proportions.
   SPIDER is explicitly unavailable in this figure: a fresh call to the pinned
   native `simulate_10X_3d` at 600,000 cells reproduced
   `ValueError: cannot reshape array of size 592704 into shape (84,84)` in
   `enhance.py:70`. This is a limitation of this implementation/run, not a
   general claim about SPIDER. The original 10,000-cell three-method overview
   and 100,000-cell scCube overview remain available in their original folders.

2. **Statistical properties** — compare expression and other statistical
   properties across methods.
   Working directory: [figure_4A](data/comparison_methods/figure_4A/).
   The historical folder name does not determine the final figure number.

3. **Cell-type clustering** —
   [clustering](data/comparison_methods/clustering/).
   **No domain clustering in Figure 5.**
   Current pipeline: all genes → normalize total to 10,000 → log1p → scale
   (clip at 10) → 30 PCs → Harmony by slice_id → 15-neighbor graph → Leiden.
   Resolution is tuned toward the known eight cell types, not maximum ARI;
   report the achieved cluster count and ARI.
   Open question from today's review: **Harmony seems not to work well here;
   is there actually a batch effect to correct?** This is a hypothesis to
   investigate, not a confirmed software failure.

4. **Spatial deconvolution** —
   [spatial_deconvolution](data/comparison_methods/spatial_deconvolution/).
   Reuse the existing ALBIS RCTD full-mode script and report correlation,
   RMSE, and the number/fraction of spots successfully scored. Spatial
   panels show slice ID 5; fitting and aggregate error metrics use all
   available input slices, subject to QC.

5. **Compute** — runtime and memory comparisons from
   [compute](data/comparison_methods/compute/).
   State dataset sizes, hardware, thread/GPU settings, and which stages
   are included (expression generation/training, spatial layout, aggregation).
   Keep simulator compute separate from downstream clustering/RCTD costs.

## Gene-count experiments to review

- Original competitor datasets: 2,000 genes, 10,000 cells in
  `data/comparison_methods/overview_full_35831924/`.
- The canonical ALBIS RCTD reference and query have 556 unique genes.
- A 556-gene subset experiment was submitted first (RCTD jobs SPIDER
  35973563 and scCube 35973564). It is a separate sensitivity analysis.
- **Preferred next comparison:** generate a fresh shared Splatter pool
  with 556 genes from the start, then regenerate SPIDER and retrain scCube.
  The pool completed and passed the 556-unique-gene check (job 35973622).
  Source directory: `data/comparison_methods/overview_native556_20260926/`.
- Generation followed by RCTD/plots: SPIDER 35973632; scCube 35973633.
  Deconvolution outputs: `spatial_deconvolution/native556_METHOD_JOBID/`.
- Cell clustering queued after those upstream jobs: SPIDER 35973656;
  scCube 35973657. Outputs: `clustering/cell_type/native556_METHOD_JOBID/`.
- These submissions retain 10,000 cells and the original spatial settings.
  Verify job completion and inspect logs before treating outputs as results.

## Thoughts and priorities for the next review

### Establish whether Harmony is appropriate

The current competitor generation wrapper does not inject an explicit
ALBIS-style slice-specific batch effect. A slice label alone is not evidence
of technical batch. Cell-type proportions may vary across slices because
of genuine spatial organization, and Harmony may remove part of that signal.

Run a paired **PCA → Leiden versus PCA → Harmony → Leiden** comparison on
identical cells/genes, with the same neighbors, seed and resolution-selection
rule. Inspect cell-type ARI and before/after embeddings, slice/type contingency
tables, and slice separation within cell types. Do not interpret slice mixing
alone as improvement. If evaluating batch correction itself, introduce a
controlled, documented technical batch perturbation and retain an unperturbed
baseline. Neither experiment is yet requested/submitted by these notes.

### Define what downstream recovery demonstrates

Good clustering/deconvolution indicates recoverable simulated signals; it is
not by itself evidence of greater biological realism or a better simulator.
Interpret these panels together with statistical properties. Equal gene count
does not equalize marker strength, expression depth, detected genes or cell-type
separability. The native 556-gene regeneration changes the expression model's
realization, rather than merely deleting features from the original matrix.

### Make deconvolution targets comparable

Competitor truth is contributing-cell fractions; ALBIS's stored truth is
captured-molecule fractions by type. Align these definitions before using
RMSE to rank methods. Current competitor references are cells that contributed
to the query spots, an optimistic reference setting. An independent-reference
experiment would test robustness. scCube's continuous expression is rounded
in separate RCTD input copies; this preprocessing choice must be disclosed.

### Keep the figure focused

Use a compact overview and one or two primary recovery summaries in the main
figure. Retain detailed per-type spatial panels, parameter searches, before/
after Harmony diagnostics and seed replicates as supporting figures. Preserve
native geometries and label occupancy settings rather than implying matched
capture footprints across methods. Final panel lettering/layout remains open.

## Overnight replication authorized and submitted

The ALBIS no-Harmony/no-batch diagnostic is already covered in Figure 3;
reuse it rather than submitting a duplicate control.

Two additional native-556 simulation seeds, 101 and 202, were submitted for
both competitors. Each seed regenerates the shared Splatter expression pool,
then SPIDER/scCube expression and spatial outputs (10,000 cells), followed
by the existing RCTD/plots and cell-type clustering. Analysis settings remain
fixed. Together with seed 20260921 this provides three simulation realizations;
inspect all outcomes and report variation rather than selecting a winner.

| Seed | Pool | SPIDER generation + RCTD | SPIDER clustering | scCube generation + RCTD | scCube clustering |
|---|---|---|---|---|---|
| 101 | 35973985 | 35973986 | 35973987 | 35973988 | 35973989 |
| 202 | 35973990 | 35973991 | 35973992 | 35973993 | 35973994 |

Generation folders: `data/comparison_methods/overview_native556_seed101/`
and `overview_native556_seed202/`. Downstream outputs use prefixes
`native556_seed101` and `native556_seed202`. Machine-readable submission
record: `data/comparison_methods/native556_replicate_jobs.json`.

Reference clarification: cell-versus-spot describes the input modalities.
An independent reference additionally requires that the reference cells are
not the cells used to construct the query spots. Current comparator pilots
use those contributing cells as reference. These new replicates retain that
protocol; independent-reference benchmarking has not been submitted. Merely
changing the whole Splatter seed can change type-expression signatures and
would introduce reference mismatch, so that experiment needs a defined
shared population model with independently sampled reference/query cells.
