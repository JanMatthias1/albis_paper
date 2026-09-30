"""Three-panel CPU benchmark: simulation time, setup-inclusive time, and memory."""
import argparse
import math
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
from matplotlib.transforms import Bbox
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
args = parser.parse_args()
data = pd.read_csv(args.root / 'measurements.csv')
data = data.loc[data.status.eq('ok')].sort_values('n_cells')
setup = pd.read_csv(args.root / 'sccube_setup.csv')
training_minutes = ' and '.join(f'{seconds / 60:.2f}' for seconds in setup.vae_training_seconds)
sizes = sorted(data.n_cells.unique())
large = max(sizes) > 600000
sampled_memory = 'generation_sampled_peak_rss_gib' in data and data.generation_sampled_peak_rss_gib.notna().any()
size_label = lambda n: f'{n/1e6:g}M' if n >= 1e6 else f'{n/1000:g}k'
data['total_seconds'] = data[['simulation_seconds', 'input_seconds',
                               'reference_seconds', 'training_seconds']].sum(axis=1, min_count=4)
out = args.root / 'figures/three_panel_with_setup'
out.mkdir(parents=True, exist_ok=True)
data.to_csv(out / 'plotted_measurements.csv', index=False)
# Figure 5 method palette, shared with clustering/cell_type/aggregate_final_batchzero.py.
methods = {'albis': ('ALBIS', '#8E63C7', 'o'),
           'spider': ('SPIDER', '#F2A65A', '^'),
           'sccube': ('scCube', '#43B7A5', 's')}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'legend.fontsize': 8, 'axes.linewidth': 0.8,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.1))
fig.subplots_adjust(left=0.065, right=0.985, top=0.72, bottom=0.30, wspace=0.37)
titles = [('A', 'Simulation time', 'VAE training and Splatter excluded'),
          ('B', 'Time including setup', 'Splatter + input + VAE training, where required'),
          ('C', 'CPU memory', 'Peak resident memory (RSS)')]
for ax, (letter, title, subtitle) in zip(axes, titles):
    ax.text(0.5, 1.19, title, transform=ax.transAxes,
            ha='center', fontsize=11, fontweight='bold')
    ax.text(0.5, 1.08, subtitle, transform=ax.transAxes, ha='center', fontsize=7.7, color='#555555')
    if large:
        ax.set_xscale('log')
        ax.set_xlim(min(sizes)/1000/1.15, max(sizes)/1000*1.15)
        ticks = [n for n in sizes if n in [10000, 100000, 600000] or n >= 1000000]
        ax.set_xticks([n/1000 for n in ticks], [size_label(n) for n in ticks])
        ax.tick_params(axis='x', labelrotation=40)
        ax.set_xlabel('Requested cells (log scale)')
    else:
        ax.set_xlim(0, max(sizes)/1000*1.045)
        ax.set_xticks([n/1000 for n in sizes])
        ax.set_xlabel('Requested cells (thousands)')
    ax.grid(axis='y', color='#DDDDDD', lw=0.6)
    ax.set_axisbelow(True)
for method, (label, color, marker) in methods.items():
    d = data.loc[data.method.eq(method)]
    axes[0].plot(d.n_cells/1000, d.simulation_seconds, color=color, marker=marker,
                 markersize=5, lw=1.7, label=label)
    axes[1].plot(d.n_cells/1000, d.total_seconds, color=color, marker=marker,
                 markersize=5, lw=1.7, label=label)
    # Memory remains segmented because batch history and measurement scope differ.
    groups = list(d.groupby('batch', sort=False)) if method == 'sccube' else [(None, d)]
    for i, (_, batch) in enumerate(groups):
        memory_values = batch.process_peak_rss_gib
        memory_label = 'scCube (cumulative)' if method == 'sccube' else label
        memory_marker, memory_style = marker, '--' if method == 'sccube' else '-'
        if method == 'sccube' and sampled_memory and batch.generation_sampled_peak_rss_gib.notna().all():
            memory_values = batch.generation_sampled_peak_rss_gib
            memory_label = 'scCube generation (sampled)'
            memory_marker, memory_style = 'D', ':'
        existing_labels = axes[2].get_legend_handles_labels()[1]
        axes[2].plot(batch.n_cells/1000, memory_values,
                     color=color, marker=memory_marker, markersize=5, lw=1.7,
                     linestyle=memory_style,
                     label=memory_label if memory_label not in existing_labels else None)
simulation_top = data.simulation_seconds.max()*1.4
axes[0].set(yscale='log', ylim=(min(1, data.simulation_seconds.min()*0.8), simulation_top),
            ylabel='Elapsed time (s, log scale)')
axes[0].yaxis.set_major_locator(FixedLocator([10**n for n in range(0, math.ceil(math.log10(simulation_top)))]))
axes[0].yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
total_top = max(1500, data.total_seconds.max()*1.4)
axes[1].set(yscale='log', ylim=(min(1, data.total_seconds.min()*0.8), total_top), ylabel='Elapsed time (s, log scale)')
axes[1].yaxis.set_major_locator(FixedLocator([10**n for n in range(0, math.ceil(math.log10(total_top)))]))
axes[1].yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
memory_top = data.process_peak_rss_gib.max()
if sampled_memory:
    memory_top = max(memory_top, data.generation_sampled_peak_rss_gib.max())
axes[2].set(ylim=(0, max(5.3, memory_top*1.1)), ylabel='Process peak RSS (GiB)')
axes[0].legend(frameon=False, loc='upper left')
axes[1].legend(frameon=False, loc='lower right')
axes[2].legend(frameon=False, loc='upper left')
fig.text(0.025, 0.956, 'Computational cost of native tissue simulation', fontsize=15, fontweight='bold')
fig.text(0.025, 0.887, f'PRELIMINARY  ·  CPU only, 1 thread  ·  556 genes  ·  seed 2025  ·  {size_label(min(sizes))}–{size_label(max(sizes))} cells',
         fontsize=9, color='#555555')
fig.text(0.025, 0.055,
         'B adds Splatter and input loading for SPIDER; Splatter, preprocessing and VAE training for scCube. ALBIS requires neither external setup step.\n'
         f'scCube setup is charged once per tissue size; the VAE was trained once per batch ({training_minutes} min). Setup-inclusive times are connected across batches; memory segments remain separate.\n' +
         ('C: scCube dashed = cumulative; dotted = generation RSS sampled every 100 ms (VAE and retained allocations resident). Splatter process excluded.\n'
          if sampled_memory else 'C: dashed scCube peaks include training and prior work within each batch; separate Splatter process excluded. No independent per-size scCube memory claim.\n') +
         'Times are elapsed wall times, excluding imports/startup, aggregation and export. Native outputs differ. One observation per size; no error bars.',
         fontsize=7.7, color='#555555', linespacing=1.6)
# Center headings over the complete panel, including axis labels, rather than
# over the plotting area alone. Use symmetric crops around that same center.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
panel_bounds = []
for ax in axes:
    headings = ax.texts[:2]
    for text in headings:
        text.set_visible(False)
    body = ax.get_tightbbox(renderer)
    center = (body.x0 + body.x1) / 2
    center_axes = ax.transAxes.inverted().transform((center, body.y0))[0]
    for text in headings:
        text.set_x(center_axes)
        text.set_visible(True)
    full = Bbox.union([body] + [text.get_window_extent(renderer) for text in headings])
    half_width = max(center - full.x0, full.x1 - center)
    bounds = Bbox.from_extents(center - half_width, full.y0,
                              center + half_width, full.y1)
    panel_bounds.append(bounds.transformed(fig.dpi_scale_trans.inverted()).expanded(1.045, 1.055))
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(out / f'compute_three_panel.{ext}', bbox_inches='tight', facecolor='white')
# Export each exact panel, including its centered headings, legend and axis labels.
# Keep the shared interpretation notes alongside the standalone files.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
panel_names = ['compute_A_simulation_time', 'compute_B_time_including_setup', 'compute_C_cpu_memory']
for ax, name, bounds in zip(axes, panel_names, panel_bounds):
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(out / f'{name}.{ext}', bbox_inches=bounds, facecolor='white')
(out / 'plot_style.json').write_text(json.dumps(dict(
    method_colors={method: values[1] for method, values in methods.items()},
    palette_source='sim_paper/code/comparison_methods/clustering/cell_type/aggregate_final_batchzero.py',
    panel_A_y_scale='log', panel_B_y_scale='log', panel_C_y_scale='linear',
    individual_panels=panel_names), indent=2)+'\n')
(out / 'README.md').write_text('# Compute figure and individual panels\n\n'
    'Combined: compute_three_panel. Individual: '+', '.join(panel_names)+'.\n'
    'Each is available as PNG, PDF and SVG. Method colors match the Figure 5 method-comparison palette.\n\n'
    + '\n\n'.join(t.get_text() for t in fig.texts)+'\n')
plt.close(fig)
print(out)
