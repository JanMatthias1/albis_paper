#!/usr/bin/env python3
"""Plot corrected scaling summaries; export all measurements and failures to CSV."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJECT = Path(__file__).resolve().parents[4]
# Match the existing comparison-method distribution figures.
METHODS = {'our_method': ('ALBIS', '#2c7fb8', 'o'),
           'splatter_spider': ('SPIDER + Splatter', '#d95f02', 's'),
           'splatter_sccube': ('scCube + Splatter', '#7570b3', '^')}


def collect(root):
    rows = []
    for path in sorted(root.glob('*/scaling_summary.json')):
        summary = json.loads(path.read_text())
        steps = {s['name']: s for s in summary['steps']}
        for name in METHODS:
            if name not in steps:
                continue
            step = steps[name]
            pipeline = [step] if name == 'our_method' else [steps.get('splatter_expression', {}), step]
            ok = all(s.get('manifest_status') == 'ok' and s.get('returncode') == 0 for s in pipeline)
            times = [s.get('time_v', {}) for s in pipeline]
            measured = all(all(k in t for k in ('wall_clock_seconds', 'max_rss_bytes', 'user_seconds', 'sys_seconds')) for t in times)
            row = dict(n_cells=summary['n_cells'], n_genes=summary['n_genes'],
                       technology=summary['technology'], scenario=summary['scenario'], seed=summary['seed'],
                       method=name, status='ok' if ok else step.get('manifest_status', 'failed'),
                       failure_reason=step.get('failure_reason') or '',
                       wall_seconds=sum(t['wall_clock_seconds'] for t in times) if ok and measured else None,
                       peak_memory_gib=max(t['max_rss_bytes'] for t in times) / 2**30 if ok and measured else None,
                       cpu_seconds=sum(t['user_seconds'] + t['sys_seconds'] for t in times) if ok and measured else None,
                       source=str(path.resolve()))
            if not ok and row['status'] == 'ok':
                row['status'] = 'prerequisite_failed'
            rows.append(row)
    if not rows:
        raise ValueError(f'No scaling summaries in {root}')
    if {(r['n_genes'], r['scenario']) for r in rows} != {(556, 'random_null')}:
        raise ValueError('Expected a homogeneous 556-gene random_null sweep')
    keys = [(r['technology'], r['n_cells'], r['method']) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError('Multiple replicates per point: aggregate explicitly before plotting')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=PROJECT / 'comparison_methods/figure_scaling/raw')
    parser.add_argument('--output-dir', type=Path, default=PROJECT / 'sim_paper/data/comparison_methods/compute/figures')
    args = parser.parse_args()
    rows = collect(args.input_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / 'compute_measurements.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
                         'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
                         'axes.linewidth': .8, 'savefig.dpi': 300, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(2, 3, figsize=(10, 6), sharey='row')
    for col, tech in enumerate(['cell', 'bin', 'spot']):
        subset = [r for r in rows if r['technology'] == tech]
        ticks = sorted({r['n_cells'] for r in subset})
        for row_idx, (metric, label, divisor) in enumerate([
                ('wall_seconds', 'Wall time (min; log scale)', 60),
                ('peak_memory_gib', 'Peak memory (GiB)', 1)]):
            ax = axes[row_idx, col]
            for name, (_, color, marker) in METHODS.items():
                points = sorted([r for r in subset if r['method'] == name], key=lambda r: r['n_cells'])
                ax.plot([r['n_cells'] / 1000 for r in points],
                        [r[metric] / divisor if r[metric] is not None else float('nan') for r in points],
                        color=color, marker=marker, markersize=4, linewidth=1.3)
                failed = [r['n_cells'] / 1000 for r in points if r['status'] != 'ok']
                ax.plot(failed, [1.025] * len(failed), transform=ax.get_xaxis_transform(),
                        linestyle='none', marker='x', color=color, markersize=6, clip_on=False)
            if row_idx == 0:
                ax.set_yscale('log')
                ax.set_title(tech.capitalize(), pad=20)
            else:
                ax.set_ylim(0, max(r['peak_memory_gib'] for r in rows if r['peak_memory_gib'] is not None) * 1.08)
                ax.set_xlabel('Requested cells (thousands)')
            if col == 0:
                ax.set_ylabel(label)
            ax.set_xticks([n / 1000 for n in ticks])
            ax.spines[['top', 'right']].set_visible(False)
            ax.grid(axis='y', color='#DDDDDD', linewidth=.5, alpha=.7)
    handles = [Line2D([], [], color=c, marker=m, label=label) for label, c, m in METHODS.values()]
    handles.append(Line2D([], [], color='#555555', marker='x', linestyle='none', label='Failed (position has no metric value)'))
    fig.legend(handles=handles, loc='upper center', ncol=4, frameon=False)
    fig.text(.5, .025, '556 genes · random_null · CPU only · one run per point\n'
             'Splatter included: times summed, peak memory = maximum across sequential stages. Lines are guides only.',
             ha='center', fontsize=8)
    fig.subplots_adjust(left=.08, right=.985, bottom=.16, top=.85, wspace=.13, hspace=.26)
    for extension in ['png', 'pdf', 'svg']:
        fig.savefig(args.output_dir / f'compute_scaling.{extension}')
    plt.close(fig)
    print(f'Saved plots and {len(rows)} measurement rows to {args.output_dir}')


if __name__ == '__main__':
    main()
