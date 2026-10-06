#!/usr/bin/env python
"""Intact 3D tissue, section-stack and single-section overview from saved data.

No simulation or expression changes. Empty/unassigned captures are hidden,
matching plot.py. PNG/PDF/SVG and a plotting summary are saved next to the inputs.
"""
import argparse
import json
from pathlib import Path
import sys
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgb
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import CELLTYPE_COLORS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--slice-number', type=int, default=5, help='One-based section number')
    parser.add_argument('--slice-id', type=int, help='Zero-based stored ID; overrides --slice-number')
    parser.add_argument('--section-spacing', type=float, default=1.25,
                        help='Display-only section spacing multiplier, matching the v12 stack')
    parser.add_argument('--output-prefix', default='overview_combined')
    args = parser.parse_args()
    if args.slice_id is not None:
        args.slice_number = args.slice_id + 1
    cfg = json.loads((args.input / 'settings.json').read_text())
    if not 1 <= args.slice_number <= cfg['n_slices']:
        parser.error('--slice-number must be between 1 and n_slices')
    if not np.isfinite(args.section_spacing) or args.section_spacing < 1:
        parser.error('--section-spacing must be finite and at least 1')
    colors = {f'type{i}': CELLTYPE_COLORS[f'Cell Type {i}'] for i in range(1, 9)}
    colors.update(cfg.get('celltype_colors', {}))
    methods = ['albis', 'sccube', 'spider']
    names = ['ALBIS', 'scCube', 'Spider']
    modalities = ['cell', 'bin', 'spot']
    datasets = {}
    for method in methods:
        if method in cfg.get('unavailable_methods', {}):
            continue
        for modality in modalities:
            a = ad.read_h5ad(args.input / method / f'{modality}.h5ad', backed='r')
            datasets[method, modality] = dict(xyz=np.asarray(a.obsm['spatial_3d']).copy(),
                sections=a.obs.slice_id.to_numpy(dtype=int),
                labels=a.obs.cell_type_true.astype(str).to_numpy(),
                empty=a.obs.is_empty.to_numpy(dtype=bool) if 'is_empty' in a.obs else np.zeros(a.n_obs, bool),
                n_genes=a.n_vars)
            if modality == 'cell' and method == 'sccube':
                native = np.asarray(a.obsm['spatial_3d_native'])
                factor = cfg['extent_um'] / cfg['sccube_grid_size']
                if not np.allclose(datasets[method, modality]['xyz'], native * factor):
                    raise ValueError('scCube coordinates differ from scaled native 3D output')
            if modality == 'cell' and len(np.unique(datasets[method, modality]['xyz'][:, 2])) <= cfg['n_slices']:
                raise ValueError(f'{method}: cell Z coordinates are section planes, not intact tissue')
            a.file.close()
            datasets[method, modality]['xyz'] += np.asarray(cfg.get('display_offsets_um', {}).get(method, [0, 0, 0]))
    lo = min(0., min(d['xyz'].min() for d in datasets.values()))
    hi = max(cfg['extent_um'], max(d['xyz'].max() for d in datasets.values()))
    if 'display_limits_um' in cfg:
        lo, hi = cfg['display_limits_um']
    pad = .1 * (hi - lo)
    marker_sizes = cfg.get('plot_marker_sizes', {'cell': 1.5, 'bin': 1.5, 'spot': 9})
    spacing = cfg['extent_um'] / cfg['n_slices'] * args.section_spacing
    planes = (np.arange(cfg['n_slices']) + .5) * spacing
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(19, 6.8), facecolor='white')
    # Paired groups each have three modalities; one row per method.
    x_tissue = .105
    x_left = [.285, .415, .545]
    x_right = [.675, .805, .935]
    centers = x_left + x_right
    y_centers = [.73, .48, .23]
    fig.text(x_tissue, .963, 'Tissue Simulation', ha='center', va='center', fontsize=15, fontweight='bold')
    fig.text(.415, .963, 'Stacked Tissue Slices', ha='center', va='center', fontsize=15, fontweight='bold')
    fig.text(x_tissue, .905, 'Intact 3D cells', ha='center', va='center', fontsize=13)
    section_title = f'Slice ID {args.slice_id}' if args.slice_id is not None else f'Slice {args.slice_number}'
    fig.text(.805, .963, section_title, ha='center', va='center', fontsize=15, fontweight='bold')
    for x, label in zip(centers, modalities * 2):
        fig.text(x, .905, label.capitalize(), ha='center', va='center', fontsize=13)
    for divider in [.21, .61]:
        fig.add_artist(plt.Line2D([divider, divider], [.115, .94], transform=fig.transFigure, color='#DDDDDD', linewidth=.8))
    summary = []
    default_limits = (lo, hi)
    for row, (method, name) in enumerate(zip(methods, names)):
        y = y_centers[row]
        name = cfg.get('method_labels', {}).get(method, name)
        fig.text(.012, y, name, rotation=90, ha='center', va='center', fontsize=12)
        if method in cfg.get('unavailable_methods', {}):
            for x in (x_tissue, .415, .805):
                fig.text(x, y, cfg['unavailable_methods'][method], ha='center', va='center', fontsize=11, color='#666666')
            continue
        lo, hi = cfg.get('method_display_limits_um', {}).get(method, default_limits)
        pad = .1 * (hi - lo)
        extent = cfg.get('method_extent_um', {}).get(method, cfg['extent_um'])
        spacing = extent / cfg['n_slices'] * args.section_spacing
        planes = (np.arange(cfg['n_slices']) + .5) * spacing
        marker_sizes = cfg.get('method_marker_sizes', {}).get(method, cfg.get(
            'plot_marker_sizes', {'cell': 1.5, 'bin': 1.5, 'spot': 9}))
        # Preserve native continuous 3D cell positions; never reconstruct from slices.
        tissue = datasets[method, 'cell']
        keep_tissue = (tissue['labels'] != 'unassigned') & ~tissue['empty']
        points = tissue['xyz'][keep_tissue]
        tissue_labels = tissue['labels'][keep_tissue]
        order = np.random.default_rng(cfg['seed']).permutation(len(points))
        ax = fig.add_axes((x_tissue-.065, y-.12, .13, .24), projection='3d')
        ax.scatter(*points[order].T, c=[colors[t] for t in tissue_labels[order]],
                   s=marker_sizes['cell'], alpha=cfg.get('plot_alpha', {}).get('cell', .9),
                   linewidths=0, depthshade=False, rasterized=True)
        ax.set_proj_type('ortho'); ax.view_init(elev=12, azim=-60)
        ax.set(xlim=(lo-pad,hi+pad), ylim=(lo-pad,hi+pad), zlim=(lo-pad,hi+pad))
        ax.set_box_aspect((1,1,1), zoom=1.55); ax.set_axis_off()
        summary.append(dict(method=method, modality='intact_tissue',
                            total_observations=len(tissue['xyz']), displayed=int(keep_tissue.sum()),
                            coordinates='Native continuous 3D cell coordinates; scCube uniform unit scaling and ALBIS display translation only',
                            section_projection=False))
        for col, modality in enumerate(modalities):
            alpha = cfg.get('plot_alpha', {}).get(modality, .75)
            section_size = cfg.get('section_marker_sizes', {}).get(modality, marker_sizes[modality])
            d = datasets[method, modality]
            xyz, sections, labels = d['xyz'], d['sections'], d['labels']
            if not np.isin(sections, np.arange(cfg['n_slices'])).all():
                raise ValueError(f'Invalid section IDs: {method}/{modality}')
            visible = (labels != 'unassigned') & ~d['empty']
            stack = xyz[visible].copy()
            stack[:, 2] = planes[sections[visible]]
            lab = labels[visible]
            order = np.random.default_rng(cfg['seed']).permutation(len(stack))
            ax = fig.add_axes((x_left[col]-.062, y-.12, .124, .24), projection='3d')
            ax.scatter(*stack[order].T, c=[colors[t] for t in lab[order]],
                       s=marker_sizes[modality], marker='s' if modality=='bin' else 'o',
                       alpha=alpha, linewidths=0, depthshade=False, rasterized=True, clip_on=False)
            ax.set_proj_type('ortho')
            ax.view_init(elev=12, azim=-60)
            ax.set(xlim=(lo-pad, hi+pad), ylim=(lo-pad, hi+pad), zlim=(0, cfg['n_slices']*spacing))
            ax.set_box_aspect((hi-lo+2*pad, hi-lo+2*pad, .65*cfg['n_slices']*spacing), zoom=1.45)
            ax.set_axis_off()

            section_mask = sections == args.slice_number - 1
            keep = visible & section_mask
            xy, lab = xyz[keep, :2], labels[keep]
            ax = fig.add_axes((x_right[col]-.06, y-.12, .12, .24))
            if modality == 'bin':
                all_xy = xyz[section_mask, :2]
                pitch = np.array([np.median(np.diff(np.unique(all_xy[:, j]))) for j in range(2)])
                if not np.all(np.isfinite(pitch) & (pitch > 0)):
                    raise ValueError(f'Cannot determine grid pitch: {method}')
                corners = np.array([[-.5,-.5],[.5,-.5],[.5,.5],[-.5,.5]])
                vertices = xy[:, None, :] + corners[None, :, :] * pitch
                faces = [alpha*np.array(to_rgb(colors[t]))+(1-alpha) for t in lab]
                ax.add_collection(PolyCollection(vertices, facecolors=faces, edgecolors='white',
                                                  linewidths=cfg.get('plot_bin_edge_width', .18), antialiaseds=False))
            else:
                order = np.random.default_rng(cfg['seed']).permutation(len(xy))
                ax.scatter(*xy[order].T, c=[colors[t] for t in lab[order]],
                           s=section_size, alpha=alpha, linewidths=0, rasterized=True)
            ax.set(xlim=(lo-pad, hi+pad), ylim=(lo-pad, hi+pad), aspect='equal')
            ax.set_axis_off()
            summary.append(dict(method=method, modality=modality, n_genes=d['n_genes'],
                                total_observations=len(xyz), stack_displayed=int(visible.sum()),
                                section_observations=int(section_mask.sum()), section_displayed=int(keep.sum())))
    handles = [plt.Line2D([], [], marker='o', ls='', markerfacecolor=c,
                         markeredgecolor='none', label=f'Cell Type {i}')
               for i,c in enumerate(colors.values(),1)]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5,.035 if cfg.get('overview_note') else .018), ncol=8,
               frameon=False, fontsize=10, handletextpad=.35, columnspacing=1.3)
    if cfg.get('overview_note'):
        fig.text(.5, .012, cfg['overview_note'], ha='center', fontsize=8)
    for ext in ('png','pdf','svg'):
        path=args.input/f'{args.output_prefix}.{ext}'
        fig.savefig(path,dpi=750,facecolor='white')
        print(path,flush=True)
    plt.close(fig)
    (args.input/f'{args.output_prefix}_summary.json').write_text(json.dumps(dict(
        input=str(args.input.resolve()), slice_number=args.slice_number, slice_id=args.slice_number-1,
        section_spacing=args.section_spacing, settings=cfg, panels=summary),indent=2)+'\n')


if __name__ == '__main__':
    main()
