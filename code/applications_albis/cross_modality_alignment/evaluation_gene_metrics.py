#!/usr/bin/env python3
"""All-gene, cross-modality grid similarity and ridge plots for Figure 4C.

Uses the benchmark's compute_grid_similarity unchanged. Recover all shared genes
from raw QC .X, restricted/reordered to the observations retained by STAIR;
normalize each observation to 10,000 counts, then log1p once. Do not use the
pre-batch layer or STAIR's filtered/normalized expression matrix.

Score bin16um/spot and bin16um/cell at unaligned, STAIR fine and truth
coordinates. A single bin16um-only rigid correction puts each stage in the truth
frame before constructing the axis-aligned grid. There is no per-target fitting.
Use the benchmark's sim cross-tech grid edges: 200 um (bin/spot), 300 um (pairs
with cell), explicitly overridden in the helper because these coordinates are
micrometers, not image pixels. Overlap grids are recomputed per pair/stage.

The _all scores use all density-qualified overlap grids; _nz uses grids with
positive expression in both modalities. All shared genes are attempted. Constant
or unsupported scores stay NaN and are omitted from densities with valid counts
shown. No overlap is undefined, never zero. Ground truth is an empirical spatial
baseline, not a theoretical upper bound for noisy cross-platform expression.

The 2x2 ridge plot shows equal-weight per-gene averages of the two reference
pairs, requiring both scores to be valid and both pairs to pass grid support.

Run via run_evaluation_gene_metrics.sh; --self-test runs small analytic checks.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import gaussian_kde

CODE = Path(__file__).resolve().parent
DEFAULT_ROOT = CODE.parents[2] / 'data/figure_4/cross_modality_alignment/independent_offsets_shift3x'
DEFAULT_HELPER = Path('/dcs04/hicks/data/multi-sample-alignment-benchmark/code/Jan/evaluation_metrics/feature_similarity.py')
TECHS = ('bin16um', 'spot', 'cell')
EVALUATION_PAIRS = (('bin16um', 'spot'), ('bin16um', 'cell'))
AVERAGE_DEFINITION = 'Per-gene arithmetic mean of bin16um–spot and bin16um–cell; equal weights; both scores required'
TECH_NAMES = {'bin16um': 'Bin (16 µm)', 'spot': 'Spot', 'cell': 'Cell'}
HELPER_TECH = {'bin16um': 'bins', 'spot': 'spots', 'cell': 'cell_obs_sectioned'}
STAGES = {'unaligned': 'spatial', 'stair_fine': 'transform_fine', 'truth': 'spatial_true'}
STAGE_LABELS = {'unaligned': 'Unaligned', 'stair_fine': 'STAIR aligned', 'truth': 'Ground truth'}
STAGE_COLORS = {'unaligned': '#B9C0C7', 'stair_fine': '#4C8FD5', 'truth': '#4A4A4A'}
METRICS = {'pcc': 'Pearson correlation', 'cos_sim': 'Cosine similarity',
           'ssim': '1-D SSIM', 'mi': 'Mutual information (nats)'}
SCORE_COLUMNS = [f'{m}_{suffix}' for suffix in ('all', 'nz') for m in METRICS]


def load_helper(path):
    spec = importlib.util.spec_from_file_location('crossmod_feature_similarity', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalize_counts(x):
    x = sparse.csr_matrix(x, dtype=np.float64)
    if not np.isfinite(x.data).all() or np.any(x.data < 0) or not np.allclose(x.data, np.round(x.data)):
        raise ValueError('QC .X must contain finite nonnegative raw counts')
    totals = np.asarray(x.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise ValueError('Retained observations must have positive library sizes')
    x = x.multiply((1e4 / totals)[:, None]).tocsr()
    x.data = np.log1p(x.data)
    return x


def fingerprint(path):
    path = Path(path).resolve()
    stat = path.stat()
    return {'path': str(path), 'size_bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def load_inputs(root, slice_id):
    from reference_metrics import rigid_fit
    aligned_path = root / f'STAIR/cross_tech/slice_{slice_id}/adata_results/Sim_CrossTech_STAIR_slice_{slice_id}.h5ad'
    aligned = ad.read_h5ad(aligned_path, backed='r')
    try:
        obs = aligned.obs.copy()
        coords = {stage: np.asarray(aligned.obsm[key], dtype=float)[:, :2]
                  for stage, key in STAGES.items()}
        aligned_genes = aligned.var_names.astype(str).tolist()
    finally:
        aligned.file.close()
    if not obs.index.is_unique or set(obs.technology.astype(str)) != set(TECHS):
        raise ValueError('Expected unique observations in all three modalities')
    for values in coords.values():
        if values.shape != (len(obs), 2) or not np.isfinite(values).all():
            raise ValueError('Invalid spatial coordinates')
    mask_ref = obs.technology.astype(str).to_numpy() == 'bin16um'
    transforms = {}
    for stage in ('unaligned', 'stair_fine'):
        rotation, shift = rigid_fit(coords[stage][mask_ref], coords['truth'][mask_ref])
        coords[stage] = coords[stage] @ rotation + shift
        transforms[stage] = {'rotation': rotation.tolist(), 'translation_um': shift.tolist()}

    datasets, files = {}, [fingerprint(aligned_path)]
    for tech in TECHS:
        mod = 'bin' if tech == 'bin16um' else tech
        path = root / 'data' / tech / f'simulation_{mod}_z_qc.h5ad'
        raw = ad.read_h5ad(path, backed='r')
        try:
            if not raw.obs_names.is_unique or not raw.var_names.is_unique:
                raise ValueError(f'{tech}: duplicate observation or gene IDs')
            mask = obs.technology.astype(str).to_numpy() == tech
            names = obs.index[mask].astype(str)
            suffix = f'_{tech}'
            if not all(name.endswith(suffix) for name in names):
                raise ValueError(f'{tech}: unexpected STAIR observation naming')
            source_names = [name[:-len(suffix)] for name in names]
            indices = raw.obs_names.get_indexer(source_names)
            if (indices < 0).any():
                raise ValueError(f'{tech}: could not match every retained observation')
            subset = raw[indices].to_memory()
        finally:
            raw.file.close()
        datasets[tech] = ad.AnnData(X=subset.X.copy(), obs=obs.loc[names].copy(), var=subset.var.copy())
        for stage in STAGES:
            datasets[tech].obsm[stage] = coords[stage][mask].copy()
        files.append(fingerprint(path))
        del subset
    genes = sorted(set.intersection(*(set(a.var_names) for a in datasets.values())))
    if not genes:
        raise ValueError('No genes shared across the three QC inputs')
    for tech, a in datasets.items():
        # Normalize against the full QC library, before selecting shared genes.
        a.X = normalize_counts(a.X)
        datasets[tech] = a[:, genes].copy()
    info = {'inputs': files, 'reference_modality': 'bin16um', 'reference_transforms': transforms,
            'n_genes_shared': len(genes), 'n_genes_stair_output': len(aligned_genes),
            'genes_recovered_from_qc': sorted(set(genes) - set(aligned_genes)),
            'n_obs': {tech: a.n_obs for tech, a in datasets.items()},
            'normalization': 'QC .X counts; per-observation total 10000 then log1p; identical across stages',
            'coordinate_frame': 'Saved STAIR preprocessed spatial_true, bin16um-only rigid correction',
            'grid_support': 'pair/stage-specific density-qualified overlap; overlap-only scores, not coverage scores'}
    return datasets, genes, info


def evaluate_pair(helper, a, b, tech_a, tech_b, stage, genes, grid_um, min_fraction):
    key = tuple(sorted((HELPER_TECH[tech_a], HELPER_TECH[tech_b])))
    c1, c2 = a.obsm[stage], b.obsm[stage]
    no_overlap = np.any(np.maximum(c1.min(0), c2.min(0)) >= np.minimum(c1.max(0), c2.max(0)))
    coverage = {'n_grids_total': 0, 'n_obs_a_in_overlap': 0, 'n_obs_b_in_overlap': 0,
                'fraction_a_in_overlap': 0., 'fraction_b_in_overlap': 0.}
    status = 'no_spatial_overlap' if no_overlap else 'ok'
    if no_overlap:
        result = pd.DataFrame(index=pd.Index(genes, name='gene'))
    else:
        labels_a, labels_b, n = helper.assign_grid_labels(a, b, stage, grid_um, min_fraction)
        coverage = {'n_grids_total': int(n),
                    'n_obs_a_in_overlap': int((labels_a >= 0).sum()),
                    'n_obs_b_in_overlap': int((labels_b >= 0).sum()),
                    'fraction_a_in_overlap': float((labels_a >= 0).mean()),
                    'fraction_b_in_overlap': float((labels_b >= 0).mean())}
        result = helper.compute_grid_similarity(a, b, tech1=HELPER_TECH[tech_a],
                    tech2=HELPER_TECH[tech_b], coord_key=stage, gene_list=genes,
                    grid_size_table={key: grid_um}, min_fraction=min_fraction)
        if n < 2:
            status = 'insufficient_overlap_grids'
        result = result.reindex(genes)
        result.index.name = 'gene'
    for col in SCORE_COLUMNS:
        if col not in result:
            result[col] = np.nan
    if 'n_grids_all' not in result or status != 'ok':
        result['n_grids_all'] = coverage['n_grids_total']
        result['n_grids_nz'] = 0
    result = result.reset_index()
    result['gene_class'] = a.var.loc[genes, 'gene_class'].to_numpy() if 'gene_class' in a.var else 'unknown'
    result['pair'] = f'{tech_a}_vs_{tech_b}'
    result['stage'] = stage
    result['status'] = status
    result['grid_size_um'] = grid_um
    record = {'pair': f'{tech_a}_vs_{tech_b}', 'stage': stage, 'status': status,
              'grid_size_um': grid_um, 'n_obs_a': a.n_obs, 'n_obs_b': b.n_obs, **coverage}
    return result, record


def average_pair_scores(table, coverage, min_grids=5):
    """Average per gene, requiring both reference pairs; never pool gene rows."""
    pairs = [f'{a}_vs_{b}' for a, b in EVALUATION_PAIRS]
    selected = table[table.pair.isin(pairs)].copy()
    if selected.duplicated(['pair', 'stage', 'gene']).any():
        raise ValueError('Duplicate pair/stage/gene rows')
    frames = []
    for stage in STAGES:
        stage_rows = selected[selected.stage == stage]
        by_pair = [stage_rows[stage_rows.pair == pair].set_index('gene').sort_index()
                   for pair in pairs]
        if any(len(frame) == 0 for frame in by_pair) or not by_pair[0].index.equals(by_pair[1].index):
            raise ValueError(f'{stage}: both pairs must contain the same gene universe')
        records = [coverage[(coverage.pair == pair) & (coverage.stage == stage)] for pair in pairs]
        if any(len(record) != 1 for record in records):
            raise ValueError(f'{stage}: expected one coverage row per pair')
        eligible = [record.iloc[0].n_grids_total >= min_grids for record in records]
        averaged = pd.DataFrame(index=by_pair[0].index)
        for metric in SCORE_COLUMNS:
            values = np.column_stack([frame[metric].to_numpy(float) for frame in by_pair])
            values[:, np.logical_not(eligible)] = np.nan
            averaged[metric + '_n_pairs'] = np.isfinite(values).sum(axis=1)
            # np.mean intentionally propagates missing scores (no one-pair fallback).
            averaged[metric] = np.mean(values, axis=1)
        averaged['stage'] = stage
        averaged['status'] = ('no_spatial_overlap' if all(record.iloc[0].status == 'no_spatial_overlap' for record in records)
                             else 'insufficient_pair_support' if not all(eligible) else 'ok')
        averaged['gene_class'] = by_pair[0]['gene_class'] if 'gene_class' in by_pair[0] else 'unknown'
        frames.append(averaged.reset_index())
    return pd.concat(frames, ignore_index=True)


def ridge_plots(table, coverage, outdir, n_genes, min_grids=5):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.titlesize': 12, 'axes.labelsize': 10,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    averaged = average_pair_scores(table, coverage, min_grids)
    averaged.to_csv(outdir / 'gene_similarity_average_per_gene.csv', index=False)
    summaries = []
    for stage in STAGES:
        for metric in SCORE_COLUMNS:
            values = averaged.loc[averaged.stage == stage, metric].dropna()
            summaries.append({'stage': stage, 'metric': metric, 'n_genes_total': n_genes,
                              'n_genes_valid': len(values), 'median': values.median(), 'mean': values.mean()})
    pd.DataFrame(summaries).to_csv(outdir / 'gene_similarity_average_summary.csv', index=False)
    plot_stages = ('stair_fine', 'truth')
    for variant in ('all', 'nz'):
        fig, axes = plt.subplots(2, 2, figsize=(9, 5.6), squeeze=False)
        for ax, metric in zip(axes.flat, METRICS):
            distributions = [averaged.loc[averaged.stage == stage, f'{metric}_{variant}'].dropna().to_numpy()
                             for stage in plot_stages]
            if metric in ('pcc', 'ssim'):
                lo, hi = -1., 1.
            elif metric == 'cos_sim':
                lo, hi = 0., 1.
            else:
                values = averaged[f'{metric}_{variant}'].dropna().to_numpy()
                lo, hi = 0., max(.1, float(values.max()) * 1.08) if len(values) else 1.
            xs = np.linspace(lo, hi, 400)
            curves = []
            for values in distributions:
                if len(values) > 2 and np.std(values) > 1e-10:
                    kde = gaussian_kde(values)
                    density = kde(xs) + kde(2*lo-xs)
                    if metric != 'mi':
                        density += kde(2*hi-xs)
                    curves.append(density)
                else:
                    curves.append(None)
            peak = max([float(d.max()) for d in curves if d is not None] or [1.])
            for index, (stage, values, density) in enumerate(zip(plot_stages, distributions, curves)):
                y = len(plot_stages)-1-index
                color = STAGE_COLORS[stage]
                ax.axhline(y, color='#DDDDDD', linewidth=.6, zorder=0)
                if density is not None:
                    heights = density / peak * .65
                    ax.fill_between(xs, y, y+heights, color=color, alpha=.65)
                    ax.plot(xs, y+heights, color=color, linewidth=1)
                if len(values):
                    median = float(np.median(values))
                    ax.plot([median, median], [y, y+.23], color=color, linewidth=1.7)
                    note = f'n={len(values)}; median={median:.3f}'
                else:
                    status = averaged.loc[averaged.stage == stage, 'status'].iloc[0]
                    note = 'No spatial overlap' if status == 'no_spatial_overlap' else 'Both pair scores required'
                ax.text(.98, y+.83, note, transform=ax.get_yaxis_transform(),
                        ha='right', va='top', fontsize=10, color='#4A4A4A')
            ax.set(xlim=(lo, hi), ylim=(-.08, 1.95), yticks=[1, 0])
            ax.set_yticklabels(['Stair\naligned' if stage == 'stair_fine' else 'Ground\ntruth'
                                for stage in plot_stages], multialignment='center', va='center')
            if metric in ('pcc', 'ssim'):
                ax.set_xticks(np.linspace(-1, 1, 5))
            ax.tick_params(axis='y', length=0)
            ax.spines[['left', 'right', 'top']].set_visible(False)
            ax.set_title(METRICS[metric], pad=12)
        fig.suptitle('Average Gene Similarity', fontsize=17, fontweight='bold')
        fig.subplots_adjust(left=.10, right=.98, bottom=.09, top=.85, wspace=.35, hspace=.48)
        for ext in ('png', 'pdf', 'svg'):
            fig.savefig(outdir / f'gene_similarity_ridges_average_{variant}.{ext}', dpi=500, facecolor='white')
        plt.close(fig)


def self_test(helper):
    # Identity and pair-order invariance use many occupied grids and varying genes.
    xy = np.array([(x+.25, y+.25) for x in range(6) for y in range(6)])
    a = ad.AnnData(X=np.column_stack([np.arange(36)+1., np.tile([1., 2., 4.], 12)]),
                   var=pd.DataFrame(index=['g1', 'g2']))
    a.obsm['truth'] = xy
    b = a.copy()
    result, _ = evaluate_pair(helper, a, b, 'bin16um', 'spot', 'truth', list(a.var_names), 1., .5)
    np.testing.assert_allclose(result[['pcc_all', 'cos_sim_all', 'ssim_all']], 1., atol=1e-10)
    reverse, _ = evaluate_pair(helper, b, a, 'spot', 'bin16um', 'truth', list(a.var_names), 1., .5)
    np.testing.assert_allclose(result[SCORE_COLUMNS], reverse[SCORE_COLUMNS], equal_nan=True)
    b.obsm['truth'] = xy + 100
    absent, coverage = evaluate_pair(helper, a, b, 'bin16um', 'spot', 'truth', list(a.var_names), 1., .5)
    assert absent[SCORE_COLUMNS].isna().all().all() and len(absent) == 2
    assert coverage['status'] == 'no_spatial_overlap'
    x = normalize_counts(sparse.csr_matrix([[1, 3], [0, 2]]))
    np.testing.assert_allclose(np.expm1(x.toarray()).sum(1), 1e4)
    constant = a.copy(); constant.X[:] = 1.
    scores, _ = evaluate_pair(helper, constant, constant, 'bin16um', 'spot', 'truth', list(a.var_names), 1., .5)
    assert scores.pcc_all.isna().all()
    frames, records = [], []
    for stage in STAGES:
        for pair, value in [('bin16um_vs_spot', .2), ('bin16um_vs_cell', .6)]:
            frames.append(pd.DataFrame({'gene': ['g1', 'g2'], 'stage': stage, 'pair': pair,
                                        **{key: [value, np.nan if pair.endswith('cell') else value]
                                           for key in SCORE_COLUMNS}}))
            records.append({'stage': stage, 'pair': pair, 'n_grids_total': 10, 'status': 'ok'})
    averaged = average_pair_scores(pd.concat(frames), pd.DataFrame(records))
    np.testing.assert_allclose(averaged.loc[averaged.gene == 'g1', SCORE_COLUMNS], .4)
    assert averaged.loc[averaged.gene == 'g2', SCORE_COLUMNS].isna().all().all()
    print('PASS: identity, pair reversal, absent-overlap NaNs, normalization, constant-gene handling, '
          'equal-weight averaging and missing-pair propagation', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--root', type=Path, default=DEFAULT_ROOT)
    parser.add_argument('--outdir', type=Path, help='Default: ROOT/gene_similarity')
    parser.add_argument('--metric-module', type=Path, default=DEFAULT_HELPER)
    parser.add_argument('--slice', type=int, default=5)
    parser.add_argument('--bin-spot-grid-um', type=float, default=200.)
    parser.add_argument('--cell-grid-um', type=float, default=300.)
    parser.add_argument('--min-fraction', type=float, default=.5)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--min-ridge-grids', type=int, default=5,
                        help='Display threshold only; all raw gene scores remain in the CSV')
    parser.add_argument('--plot-only', action='store_true', help='Redraw saved CSVs without recomputing scores')
    args = parser.parse_args()
    if not all(np.isfinite(v) and v > 0 for v in
               (args.bin_spot_grid_um, args.cell_grid_um, args.min_fraction)):
        parser.error('Grid sizes and min-fraction must be finite and positive')
    if args.min_ridge_grids < 2:
        parser.error('min-ridge-grids must be at least 2')
    outdir = args.outdir or args.root / 'gene_similarity'
    if args.plot_only:
        table = pd.read_csv(outdir / 'gene_similarity_per_gene.csv')
        coverage = pd.read_csv(outdir / 'gene_similarity_coverage.csv')
        provenance = json.loads((outdir / 'provenance.json').read_text())
        ridge_plots(table, coverage, outdir, provenance['n_genes_shared'], args.min_ridge_grids)
        provenance['ridge_min_grids'] = args.min_ridge_grids
        provenance['average_definition'] = AVERAGE_DEFINITION
        provenance['plotted_pairs'] = [list(pair) for pair in EVALUATION_PAIRS]
        provenance['plot_script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        (outdir / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
        with (outdir / 'METHODS.txt').open('a') as fh:
            fh.write('\nAverage: ' + AVERAGE_DEFINITION + '\n')
            fh.write(f'\nRidge rendering requires at least {args.min_ridge_grids} overlap grids; raw CSV scores are retained.\n')
        print(f'Redrawn ridge plots in {outdir.resolve()}', flush=True)
        return
    helper = load_helper(args.metric_module)
    if args.self_test:
        self_test(helper)
        return
    outdir.mkdir(parents=True, exist_ok=True)
    datasets, genes, provenance = load_inputs(args.root, args.slice)
    print(f'Using {len(genes)} shared QC genes and retained observations: {provenance["n_obs"]}', flush=True)
    (outdir / 'gene_list.txt').write_text('\n'.join(genes) + '\n')
    frames, records = [], []
    for tech_a, tech_b in EVALUATION_PAIRS:
        grid = args.cell_grid_um if 'cell' in (tech_a, tech_b) else args.bin_spot_grid_um
        for stage in STAGES:
            print(f'{tech_a} vs {tech_b}: {stage}, grid={grid:g} µm', flush=True)
            frame, record = evaluate_pair(helper, datasets[tech_a], datasets[tech_b], tech_a, tech_b,
                                          stage, genes, grid, args.min_fraction)
            frames.append(frame); records.append(record)
            print(f'  {record["status"]}; overlap grids={record["n_grids_total"]}', flush=True)
    table, coverage = pd.concat(frames, ignore_index=True), pd.DataFrame(records)
    table.to_csv(outdir / 'gene_similarity_per_gene.csv', index=False)
    coverage.to_csv(outdir / 'gene_similarity_coverage.csv', index=False)
    summary = []
    for (pair, stage), group in table.groupby(['pair', 'stage'], sort=False):
        for metric in SCORE_COLUMNS:
            values = group[metric].dropna()
            summary.append({'pair': pair, 'stage': stage, 'metric': metric, 'n_genes_total': len(group),
                            'n_genes_valid': len(values), 'median': values.median(), 'mean': values.mean()})
    pd.DataFrame(summary).to_csv(outdir / 'gene_similarity_summary.csv', index=False)
    provenance.update({'slice_id': args.slice, 'script': str(Path(__file__).resolve()),
                       'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       'metric_module': str(args.metric_module.resolve()),
                       'metric_module_sha256': hashlib.sha256(args.metric_module.read_bytes()).hexdigest(),
                       'min_fraction': args.min_fraction, 'coverage': records,
                       'ridge_min_grids': args.min_ridge_grids,
                       'average_definition': AVERAGE_DEFINITION,
                       'plotted_pairs': [list(pair) for pair in EVALUATION_PAIRS],
                       'ridge_definition': 'gene-score KDE, Scott bandwidth, reflected at metric boundaries; '
                       'peak scaling common within each panel; median ticks; no gene selection'})
    (outdir / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    (outdir / 'METHODS.txt').write_text(__doc__ + '\n\nMetric implementation: ' + str(args.metric_module.resolve()) +
        '\nAll genes are attempted, including noise/shared markers. Valid gene counts can differ by metric/stage.\n'
        'Use gene_similarity_coverage.csv alongside the scores: the metric evaluates only overlapping grids.\n'
        'The _nz output is a sensitivity view, without an extra minimum nonzero-grid filter.\n'
        'These are gene distributions for one simulated section, not independent simulation replicates.\n'
        f'Ridge rendering requires at least {args.min_ridge_grids} overlap grids; raw CSV scores are retained.\n'
        'Average: ' + AVERAGE_DEFINITION + '\n')
    ridge_plots(table, coverage, outdir, len(genes), args.min_ridge_grids)
    print(f'Saved metrics and ridge plots to {outdir.resolve()}', flush=True)


if __name__ == '__main__':
    main()
