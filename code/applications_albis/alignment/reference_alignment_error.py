"""Evaluate saved STAIR reconstructions using a reference-only rigid frame fit.

No scale, reflection, per-target fit, or observation matching is estimated.
Observation correspondence is inherited from the saved simulation output.
"""
import argparse
import json
from pathlib import Path
import sys
import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import apply_style, BOX_EDGE, BOX_MEDIAN, GRID_COLOR
apply_style()
BASE = Path(__file__).resolve().parents[3] / 'data/figure_4/alignment'


def fit_rigid(source, target):
    """Row-vector map source @ R + t -> target; proper rotation only."""
    sx, sy = source.mean(0), target.mean(0)
    x, y = source-sx, target-sy
    if np.linalg.matrix_rank(x) < 2:
        raise ValueError('Reference coordinates do not span two dimensions')
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.eye(2)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    return rotation, sy-sx @ rotation


def evaluate(path, dataset, reference='0'):
    with h5py.File(path) as f:
        obs = read_elem(f['obs'])
        coords = {key: np.asarray(read_elem(f['obsm'][key]), dtype=float)[:, :2]
                  for key in ('spatial_true', 'spatial', 'transform_fine')}
    if not obs.index.is_unique:
        raise ValueError('Observation IDs must be unique')
    if any(a.shape != (len(obs), 2) or not np.isfinite(a).all() for a in coords.values()):
        raise ValueError('Invalid coordinate shape or nonfinite coordinates')
    ids = obs['slice_id'].astype(str).to_numpy()
    mask = ids == reference
    if mask.sum() < 3:
        raise ValueError(f'Missing/insufficient reference slice {reference}')
    truth = coords['spatial_true']
    rows, fits = [], []
    for stage, key in [('Unaligned', 'spatial'), ('STAIR aligned', 'transform_fine')]:
        # Map truth into each condition's reference frame, exactly as in the
        # supplied pairwise metric. Rigid transforms preserve distance units.
        rotation, translation = fit_rigid(truth[mask], coords[key][mask])
        expected = truth @ rotation + translation
        squared = np.sum((coords[key]-expected)**2, axis=1)
        fits.append(dict(dataset=dataset, condition=stage, reference_slice=reference,
                         rotation=rotation.tolist(), translation_um=translation.tolist()))
        for sid in sorted(set(ids), key=float):
            selected = ids == sid
            mse = float(squared[selected].mean())
            rows.append(dict(dataset=dataset, condition=stage, slice_id=sid,
                             reference_slice=reference, is_reference=sid == reference,
                             n_obs=int(selected.sum()), mse_um2=mse, rmse_um=np.sqrt(mse)))
    return rows, fits


def plot(df, out):
    targets = df[~df.is_reference]
    names = list(dict.fromkeys(targets.dataset))
    labels = {'cell_r6000':'Cell (larger sphere)', 'cell':'Cell','bin16um':'Bin (16 µm)','spot':'Spot'}
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {'Unaligned':'#B9C0C7','STAIR aligned':'#4C8FD5'}
    rng = np.random.default_rng(0)
    for i, ds in enumerate(names):
        sub = targets[targets.dataset == ds]
        jitter = dict(zip(sorted(sub.slice_id.unique(), key=float), rng.uniform(-.035,.035,sub.slice_id.nunique())))
        for stage, offset in [('Unaligned',-.18),('STAIR aligned',.18)]:
            values = sub[sub.condition == stage]
            bp = ax.boxplot(values.rmse_um, positions=[i+offset], widths=.28,
                            patch_artist=True, showfliers=False,
                            medianprops={'color':BOX_MEDIAN},
                            boxprops={'edgecolor':BOX_EDGE},
                            whiskerprops={'color':BOX_EDGE},capprops={'color':BOX_EDGE})
            bp['boxes'][0].set(facecolor=colors[stage], alpha=.65)
            ax.scatter([i+offset+jitter[s] for s in values.slice_id],values.rmse_um,
                       s=24,color=colors[stage],edgecolors=BOX_EDGE,linewidths=.5,zorder=3)
    ax.set_xticks(range(len(names)), [labels[n] for n in names])
    ax.set_ylabel('Reference-corrected alignment RMSE (µm)')
    ax.set_title('Alignment accuracy', fontsize=17, fontweight='bold', pad=14)
    # Log scale retains visibility of aligned residuals and large perturbations.
    if (targets.rmse_um <= 0).any():
        ax.set_yscale('symlog', linthresh=1)
    else:
        ax.set_yscale('log')
    ax.grid(axis='y',color=GRID_COLOR,alpha=.6)
    ax.set_axisbelow(True)
    ax.spines[['top','right']].set_visible(False)
    ax.legend(handles=[plt.Rectangle((0,0),1,1,color=c,alpha=.65,label=s) for s,c in colors.items()],frameon=False)
    fig.tight_layout()
    fig.savefig(out,dpi=300)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=Path,default=BASE)
    parser.add_argument('--reference',default='0')
    args=parser.parse_args()
    rows, fits=[],[]
    for ds in ['cell','bin16um','spot']:  # 2026-09-28: Fig4C strong-mix data, cell = 600k dense
        r,f=evaluate(args.base/f'STAIR/{ds}/adata_results/Sim_3D_STAIR_{ds}.h5ad',ds,args.reference)
        rows.extend(r);fits.extend(f)
    out=args.base/'plots';out.mkdir(exist_ok=True)
    df=pd.DataFrame(rows)
    df.to_csv(out/'figure4a_reference_alignment_error.csv',index=False)
    (out/'figure4a_reference_alignment_transforms.json').write_text(json.dumps(fits,indent=2))
    plot(df,out/'figure4a_reference_alignment_error.png')
    print(df[~df.is_reference].groupby(['dataset','condition']).rmse_um.agg(['count','median','min','max']).to_string())
    print('Reference fit diagnostics:')
    print(df[df.is_reference][['dataset','condition','rmse_um']].to_string(index=False))


if __name__=='__main__':
    main()
