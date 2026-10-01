# Current compute figures and native benchmark

## Three-seed run (2026-10-01, albis 0.1.2)

`bash submit_seeds.sh` prepares `data/figure_5/compute/native_seeds_<date>/` and submits
60 jobs: per seed (2025, 101, 202) one Splatter reference, ALBIS and SPIDER at each of
10k/100k/200k/400k/600k/1M/2M/5M cells (`native_extension.sbatch`), and three scCube
batches (10k-100k, 200k-600k, 1M-5M, one VAE each; `sccube_seed_batch.sbatch`, 100 ms
RSS sampling). `prepare_seeds.py` writes every size x seed from
`settings_template_n600000_seed2025.json` with the same (n/600000)^(1/3) scaling as
before (its seed-2025/101/202 files match the earlier pilot's exactly) and freezes a
code snapshot; jobs run that snapshot. `sccube_split.py` gained `--seed` and `--out`
(defaults unchanged). Job IDs: `<run>/submissions.tsv`. Current run:
`native_seeds_20261001`. The report/plot below still reads the single-seed run's
table; a multi-seed version (mean and spread over seeds) is the next step.
The single-seed run (albis 0.1.1) is archived in
`data/figure_5/archive_albis0.1.1_20261001/compute/`.

This directory retains the code for the current preliminary four-panel figure,
its four standalone images, and the completed native scaling workflow through
5 million cells. Additional runs and seed replicates are deferred pending feedback.

## Reproduce the current images

From `/dcs04/hicks/data/Jan/sim_project`:

```bash
MPLCONFIGDIR=/tmp/compute-preview-mpl OPENBLAS_NUM_THREADS=1 \
comparison_methods/env/analysis/bin/python \
sim_paper/code/comparison_methods/compute/plot_compute_four_panel_preliminary.py \
--root sim_paper/data/figure_5/compute/native_large_1m_5m_20260929
```

This reads the existing 24-measurement table; it does not launch simulations or
change the underlying measurements. Outputs are under the run's
`figures/four_panel_preliminary/`, each in PNG, PDF and SVG:

- `compute_four_panel`: combined figure.
- `simulation_time`: generation time excluding setup.
- `time_including_setup`: time with setup charged at each size.
- `memory_during_generation`: available generation measurements and marked proxies.
- `memory_including_setup`: partial method-process peaks; Splatter RAM unmeasured.

Colors match Figure 5. Titles and subtitles are centered on each complete image;
individual exports exclude the shared figure caption. The plot folder retains
interpretation notes and plotted measurements. `METHODS.md` provides manuscript
text, also saved beside the images.

## Retained run code

- `albis_native.py`, `spider_native.py`, `sccube_split.py`: direct native calls,
  timing and read-only validation; scCube reuses an in-memory VAE within each batch.
- `reference_native.R`, `reference_native.sbatch`: Splatter reference generation.
- `native_extension.sbatch`: independent ALBIS/SPIDER jobs.
- `sccube_large.sbatch`, `sample_rss.py`: scCube large-size batch and 100 ms RSS sampling.
- `prepare_large.py`: prepare 1M/2M/5M inputs and freeze source snapshots.
- `report_large.py`, `report_large.sbatch`: combine prior and large-size results.
  The canonical report now invokes the current preliminary four-panel plot.

The current completed run is
`sim_paper/data/figure_5/compute/native_large_1m_5m_20260929`.
It includes eight cell counts (10k, 100k, 200k, 400k, 600k, 1M, 2M, 5M),
three methods, 556 genes, one computational CPU thread and seed 2025.
Earlier source runs remain under the data compute directory's
`_archive/superseded_20260929/`. The active preparation/report scripts resolve
those archived inputs. Existing run `code/` directories are immutable historical
snapshots; the retained Slurm launchers execute those snapshots. Use the canonical
plot command above for the current figure styling.

## Interpretation

Native outputs differ: ALBIS produces explicit molecules, scCube produces
expression and metadata, and SPIDER retains its native reference-backed view.
No output transformations or model checkpoints are introduced.

Time including setup charges the batch's measured setup cost at every size;
scCube was not independently retrained for each point. Generation-memory samples
exist for scCube only at 1M–5M. ALBIS/SPIDER whole-process peaks are explicit
proxies. The setup-memory panel excludes unmeasured Splatter-process memory and
includes earlier within-batch work for scCube. These are preliminary, unmatched
memory scopes, not fully isolated with/without-setup measurements.

## Archive

`_archive/superseded_20260929/` preserves older figures, smoke/pilot/extension
entry points, legacy plans and the unsubmitted fresh-process memory workflow.
Its manifest records original locations and hashes. Nothing was deleted;
existing raw results, detailed timing tables and executed source snapshots remain
unchanged. The fresh-process workflow and additional seeds are not active.
