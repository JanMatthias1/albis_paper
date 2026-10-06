"""Four-panel compute figure from report_seeds.py's measurements.csv: lines = mean over
simulation seeds, small points = individual seeds. Writes <root>/figures/compute/."""
import argparse
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/compute-figure-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
args = parser.parse_args()
root = args.root.resolve()
out = root / 'figures' / 'compute'
out.mkdir(parents=True, exist_ok=True)

data = pd.read_csv(root / 'measurements.csv')
failed = data.loc[data.status.eq('failed')]
data = data.loc[data.status.eq('ok')]
n_seeds = data.seed.nunique()
methods = {'albis': ('ALBIS', '#8E63C7', 'o'), 'sccube': ('scCube', '#43B7A5', 's'), 'spider': ('SPIDER', '#F2A65A', '^')}
panels = [
    ('simulation_seconds', 'Simulation time', '3D tissue to 10 cell-level sections; Splatter and scCube VAE training excluded',
     'Elapsed time (s)', 'simulation_time'),
    ('total_seconds', 'Time including setup', 'Splatter reference + input loading + scCube VAE training',
     'Elapsed time (s)', 'time_including_setup'),
    ('generation_peak_rss_gib', 'Memory during generation', 'ALBIS, SPIDER: one size per process; scCube: 100 ms samples',
     'Peak RAM (GiB)', 'memory_during_generation'),
    ('setup_peak_rss_gib', 'Memory including setup', 'Larger of method-process and Splatter-process peaks',
     'Peak RAM (GiB)', 'memory_including_setup')]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 500, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
fig.subplots_adjust(left=.10, right=.97, bottom=.19, top=.81, wspace=.33, hspace=.85)
sizes = sorted(data.n_cells.unique())
ticks, labels = [10000, 100000, 600000, 1000000, 2000000, 5000000], ['10k', '100k', '600k', '1M', '2M', '5M']
for ax, (key, title, subtitle, ylabel, name) in zip(axes.flat, panels):
    ax.text(.5, 1.22, title, transform=ax.transAxes, ha='center', fontsize=12, fontweight='bold')
    ax.text(.5, 1.10, subtitle, transform=ax.transAxes, ha='center', fontsize=8, color='#555555')
    for method, (label, color, marker) in methods.items():
        d = data.loc[data.method.eq(method)]
        if key == 'setup_peak_rss_gib' and method == 'sccube':
            label += ' (includes training)'
        mean = d.groupby('n_cells')[key].mean().reindex(sizes)  # reindex shows gaps for missing sizes
        ax.scatter(d.n_cells, d[key], s=10, color=color, alpha=.35, linewidths=0, zorder=2)
        ax.plot(sizes, mean, label=label, color=color, marker=marker, lw=1.7, markersize=5, zorder=3,
                markerfacecolor=color, linestyle='-')
    ax.set_xscale('log')
    ax.set_xticks(ticks, labels, rotation=35)
    ax.set_xlim(8500, 6000000)
    ax.set_xlabel('Requested cells (log scale)')
    ax.set_ylabel(ylabel)
    if key.endswith('seconds'):
        ax.set_yscale('log')
    else:
        ax.set_ylim(bottom=0)
    ax.grid(axis='y', color='#DDDDDD', lw=.6)
    ax.legend(frameon=False, fontsize=8)
fig.text(.5, .975, 'Computational cost of native tissue simulation', ha='center', fontsize=15, fontweight='bold')
fig.text(.5, .945, f'Mean of {n_seeds} simulation seeds (points: individual seeds) · CPU, 1 thread · 556 genes',
         ha='center', color='#555555')
notes = ('Every method runs each size in its own process; scCube trains its own VAE in each run (setup is measured per run).\n'
         'Memory is the process peak RAM; scCube generation-only memory is the peak of 100 ms samples during generation.\n'
         'Native output representations differ: ALBIS generates explicit molecules; SPIDER and scCube generate cell-level expression.')
# Sizes where a method failed for every seed: no point is drawn; say so under the figure.
for (method, n), runs in failed.groupby(['method', 'n_cells']):
    if not data.loc[data.method.eq(method) & data.n_cells.eq(n)].empty:
        continue
    stage = runs.failed_stage.dropna().unique() if 'failed_stage' in runs else []
    notes += (f"\n{methods[method][0]}, {n:,} cells: did not complete ({len(runs)} of {len(runs)} seeds failed"
              + (f" in {', '.join(stage)}" if len(stage) else '') + '; see RESULTS.md).')
fig.text(.07, .025, notes, fontsize=8, linespacing=1.5)

# Save each panel on its own as well, with its headings centred on the panel body.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
bounds = []
for ax in axes.flat:
    headings = ax.texts[:2]
    for t in headings:
        t.set_visible(False)
    body = ax.get_tightbbox(renderer)
    center = (body.x0 + body.x1) / 2
    for t in headings:
        t.set_x(ax.transAxes.inverted().transform((center, body.y0))[0])
        t.set_visible(True)
    full = Bbox.union([body] + [t.get_window_extent(renderer) for t in headings])
    width = max(center - full.x0, full.x1 - center)
    bounds.append(Bbox.from_extents(center - width, full.y0, center + width, full.y1)
                  .transformed(fig.dpi_scale_trans.inverted()).expanded(1.045, 1.055))
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(out / f'compute_four_panel.{ext}', bbox_inches='tight', facecolor='white')
    for text in fig.texts:
        text.set_visible(False)
    # Hide the other panels too: bbox_inches only crops the view, so PDF/SVG would still
    # contain all four panels (they show up when the file is imported into an editor).
    for i, (panel, bound) in enumerate(zip(panels, bounds)):
        for j, other in enumerate(axes.flat):
            other.set_visible(i == j)
        fig.savefig(out / f'{panel[-1]}.{ext}', bbox_inches=bound, facecolor='white')
    for ax in axes.flat:
        ax.set_visible(True)
    for text in fig.texts:
        text.set_visible(True)
plt.close(fig)
print(f'{len(data)} measurements from {n_seeds} seeds; figures in {out}', flush=True)
