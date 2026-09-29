"""Manuscript layout previews from recorded measurements, with a conceptual output schematic."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
args = parser.parse_args()
root = args.root
data = pd.read_csv(root / 'measurements.csv')
data = data.loc[data.status.eq('ok')].sort_values('n_cells')
setup = pd.read_csv(root / 'sccube_setup.csv')
out = root / 'figures/panel_preview'
out.mkdir(parents=True, exist_ok=True)
colors = {'albis': '#4C8FD5', 'spider': '#8E63C7', 'sccube': '#8A949E'}
labels = {'albis': 'ALBIS', 'spider': 'SPIDER', 'sccube': 'scCube'}
markers = {'albis': 'o', 'spider': '^', 'sccube': 's'}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.titlesize': 11, 'axes.labelsize': 9,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'legend.fontsize': 8, 'axes.linewidth': 0.8,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300, 'svg.fonttype': 'none', 'pdf.fonttype': 42})


def heading(ax, letter, title, subtitle):
    ax.text(-0.13, 1.16, letter, transform=ax.transAxes, fontsize=17, fontweight='bold')
    ax.text(0, 1.16, title, transform=ax.transAxes, fontsize=11, fontweight='bold')
    ax.text(0, 1.07, subtitle, transform=ax.transAxes, fontsize=8.2, color='#555555')


def tissue_axis(ax):
    ax.set_xlim(0, 625)
    ax.set_xticks([10, 100, 200, 400, 600])
    ax.set_xlabel('Requested cells (thousands)')
    ax.grid(axis='y', color='#DDDDDD', linewidth=0.6)
    ax.set_axisbelow(True)


def runtime(ax):
    heading(ax, 'A', 'Simulation runtime', 'Prepared reference / VAE; setup excluded')
    for method in ['albis', 'spider', 'sccube']:
        d = data.loc[data.method.eq(method)]
        ax.plot(d.n_cells / 1000, d.simulation_seconds, color=colors[method],
                marker=markers[method], markersize=5, linewidth=1.7, label=labels[method])
    tissue_axis(ax)
    ax.set_ylim(0, 125)
    ax.set_ylabel('Simulation time (s)')
    ax.legend(frameon=False, loc='upper left')


def setup_panel(ax):
    heading(ax, 'B', 'Setup measured separately', 'Individual measurements; no tissue-size axis')
    rows = [
        ('Splatter reference', [float(setup.reference_generation_seconds.iloc[0])], '#B9C0C7'),
        ('SPIDER input / preprocessing', data.loc[data.method.eq('spider'), 'input_seconds'].tolist(), colors['spider']),
        ('scCube input / preprocessing', setup.input_seconds.tolist(), colors['sccube']),
        ('scCube VAE training', setup.vae_training_seconds.tolist(), colors['sccube']),
    ]
    for i, (name, values, color) in enumerate(rows):
        offsets = np.linspace(-0.13, 0.13, len(values)) if len(values) > 1 else [0.]
        for value, delta in zip(values, offsets):
            ax.scatter(value, i + delta, color=color, s=35, edgecolors='white', linewidth=0.5, zorder=3)
        if i == 0:
            ax.annotate(f'{values[0]:.2f} s', (values[0], i), xytext=(7, 0),
                        textcoords='offset points', va='center', fontsize=8)
        if i == 3:
            for value, delta in zip(values, offsets):
                ax.annotate(f'{value/60:.2f} min', (value, i + delta),
                            xytext=(-7, -8 if delta > 0 else 8), textcoords='offset points',
                            va='center', ha='right', fontsize=8)
    ax.set_xscale('log')
    ax.set_xlim(0.1, 2000)
    ax.set_ylim(3.65, -0.6)
    ax.set_yticks(range(4), [r[0] for r in rows])
    ax.xaxis.set_major_locator(FixedLocator([0.1, 1, 10, 100, 1000]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    ax.set_xlabel('Setup time (s, log scale)')
    ax.grid(axis='x', color='#DDDDDD', linewidth=0.6)
    ax.set_axisbelow(True)
    ax.text(0, -0.27, 'ALBIS: no external reference or VAE training.\nVAE: two training runs on different nodes; one training per batch.',
            transform=ax.transAxes, fontsize=7.5, color='#555555', va='top')


def memory(ax):
    heading(ax, 'C', 'Memory footprint', 'Solid: independent processes; dashed: cumulative scCube RSS')
    for method in ['albis', 'spider']:
        d = data.loc[data.method.eq(method)]
        ax.plot(d.n_cells / 1000, d.process_peak_rss_gib, color=colors[method],
                marker=markers[method], markersize=5, linewidth=1.7, label=labels[method])
    for i, (_, d) in enumerate(data.loc[data.method.eq('sccube')].groupby('batch', sort=False)):
        ax.plot(d.n_cells / 1000, d.process_peak_rss_gib, 's--', color=colors['sccube'],
                markersize=5, linewidth=1.7, label='scCube (cumulative)' if i == 0 else None)
    tissue_axis(ax)
    ax.set_ylim(0, 5.3)
    ax.set_ylabel('Process peak RSS (GiB)')
    ax.legend(frameon=False, loc='upper left')
    ax.text(0, -0.25, 'scCube peaks include training and prior work within each batch.\nThe two scCube batches are deliberately not connected.',
            transform=ax.transAxes, fontsize=7.5, color='#555555', va='top')


def representation(ax):
    heading(ax, 'D', 'What each native call returns', 'Conceptual schematic; these are different output representations')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    rng = np.random.default_rng(42)
    descriptions = [
        ('albis', 0.79, 'Explicit molecules', 'Molecule XYZ + gene IDs; tissue cell metadata'),
        ('sccube', 0.47, 'Cell expression', 'Expression matrix + cell metadata / XYZ'),
        ('spider', 0.15, 'Native expression view', 'Reference-backed expression view + cell XYZ'),
    ]
    for method, y, title, detail in descriptions:
        color = colors[method]
        ax.add_patch(Rectangle((0, y-0.12), 1, 0.27, facecolor=color, alpha=0.065, edgecolor='none'))
        ax.text(0.04, y+0.08, labels[method], color=color, fontweight='bold', fontsize=10)
        if method == 'albis':
            points = rng.uniform([0.045, y-0.09], [0.245, y+0.035], size=(48, 2))
            ax.scatter(points[:, 0], points[:, 1], s=5, color=color, alpha=0.8)
        else:
            # Decorative matrix glyph only; no simulated model data are altered.
            for row in range(4):
                for col in range(7):
                    ax.add_patch(Rectangle((0.045+col*0.027, y-0.092+row*0.03), 0.023, 0.025,
                                            facecolor=color, alpha=0.2+0.1*((row+col)%6), edgecolor='none'))
            if method == 'spider':
                ax.add_patch(Rectangle((0.04, y-0.098), 0.199, 0.127,
                                        fill=False, edgecolor=color, linestyle='--', linewidth=1))
        ax.text(0.3, y+0.035, title, fontsize=10, fontweight='bold')
        ax.text(0.3, y-0.035, detail, fontsize=7.8, color='#555555')
        if method == 'albis':
            d = data.loc[data.method.eq('albis')].sort_values('n_cells')
            ax.text(0.3, y-0.09, f'{d.n_molecules.iloc[0]/1e6:.2f}M → {d.n_molecules.iloc[-1]/1e6:.2f}M molecules (10k → 600k cells)', fontsize=7.8, color=color)
    ax.text(0, -0.25, 'Cost of native representations, not equivalent outputs.\nSPIDER’s view is retained; no forced materialization is timed.',
            transform=ax.transAxes, fontsize=7.5, color='#555555', va='top')


def save(fig, name):
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(out / f'{name}.{ext}', bbox_inches='tight', facecolor='white')
    plt.close(fig)


fig, axes = plt.subplots(2, 2, figsize=(13.5, 10))
fig.subplots_adjust(left=0.085, right=0.97, bottom=0.17, top=0.85, wspace=0.76, hspace=0.95)
runtime(axes[0, 0]); setup_panel(axes[0, 1]); memory(axes[1, 0]); representation(axes[1, 1])
fig.text(0.045, 0.966, 'Computational cost of native tissue simulation', fontsize=16, fontweight='bold')
fig.text(0.045, 0.93, 'PRELIMINARY  ·  556 genes  ·  1 CPU thread  ·  seed 2025  ·  10k–600k cells', fontsize=10, color='#555555')
fig.text(0.045, 0.03,
         'One observation per method and size; no seed replicates or error bars. Lines guide the eye.\n'
         'Runtime excludes imports/startup, aggregation and export. scCube uses an in-memory VAE; no model files are saved.\n'
         'Memory excludes the separate Splatter process. scCube memory is cumulative, not an independent per-size measurement.',
         fontsize=8, color='#555555', linespacing=1.55)
save(fig, 'compute_four_panel')

# A compact alternative when the native-output explanation belongs in the caption.
fig, axes = plt.subplots(1, 3, figsize=(16, 4.9))
fig.subplots_adjust(left=0.055, right=0.985, bottom=0.3, top=0.73, wspace=0.95)
runtime(axes[0]); setup_panel(axes[1]); memory(axes[2])
fig.text(0.025, 0.945, 'Computational cost of native tissue simulation', fontsize=15, fontweight='bold')
fig.text(0.025, 0.875, 'PRELIMINARY  ·  556 genes  ·  1 CPU thread  ·  seed 2025', fontsize=9, color='#555555')
fig.text(0.025, 0.035,
         'Native outputs differ: ALBIS explicit molecules (0.88M–52.82M); scCube cell expression; SPIDER reference-backed expression view.\n'
         'One observation per size; no error bars. Runtime excludes imports/startup, aggregation and export. scCube VAE is held in memory; no models saved.\n'
         'Memory excludes separate Splatter generation; scCube memory is cumulative within each batch, not independent per size.',
         fontsize=8, color='#555555', linespacing=1.5)
save(fig, 'compute_three_panel')
print(out)
