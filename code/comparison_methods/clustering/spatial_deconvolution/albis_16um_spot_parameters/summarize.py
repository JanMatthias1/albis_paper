"""Validate ALBIS RCTD results and plot native spot truth versus estimates."""
import json
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from generate import OUT, SEEDS


def main():
    dest = OUT / 'figures'
    dest.mkdir(exist_ok=True)
    rows = []
    types = [f'type{i}' for i in range(1, 9)]
    for seed in SEEDS:
        run = OUT / f'seed{seed}'
        meta = json.loads((run / 'results/metrics_summary.json').read_text())
        pred = pd.read_csv(run / 'results/estimated_fractions_wide.csv')
        idx = pred.spot_id.str.removeprefix('spot_').astype(int).to_numpy() - 1
        a = ad.read_h5ad(run / 'inputs/spot.h5ad', backed='r')
        assert pred.spot_id.is_unique and (idx >= 0).all() and (idx < a.n_obs).all()
        truth = np.asarray(a.obsm['cell_type_frac_true'])[idx]
        est = pred[types].to_numpy()
        assert np.allclose(truth.sum(1), 1) and np.allclose(est.sum(1), 1)
        rmse = float(np.sqrt(np.mean((truth-est)**2)))
        corr = float(np.corrcoef(truth.ravel(), est.ravel())[0, 1])
        assert abs(rmse-meta['overall_rmse']) < 1e-4
        assert abs(corr-meta['overall_pearson_r']) < 1e-4
        rows.append(dict(seed=seed, rmse=rmse, pearson_r=corr, mean_absolute_error=float(np.abs(truth-est).mean()), spots_supplied=a.n_obs, spots_scored=len(idx)))
        keep = a.obs.iloc[idx].slice_id.astype(int).to_numpy() == 4
        xy = np.asarray(a.obsm['spatial'])[idx][keep, :2]
        fig, axes = plt.subplots(2, 8, figsize=(17, 5), layout='constrained')
        for j, name in enumerate(types):
            for i, values in enumerate([truth[keep], est[keep]]):
                ax = axes[i, j]
                im = ax.scatter(xy[:, 0], xy[:, 1], c=values[:, j], cmap='viridis', vmin=0, vmax=1, marker='o', s=9, linewidths=0, rasterized=True)
                ax.set_aspect('equal')
                ax.set_xticks([])
                ax.set_yticks([])
                if i == 0:
                    ax.set_title(name)
                if j == 0:
                    ax.set_ylabel('True' if i == 0 else 'RCTD')
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=.7, label='Molecule fraction')
        fig.suptitle(f'ALBIS · Slice 5 · Simulation seed {seed}')
        for ext in ['png', 'pdf', 'svg']:
            fig.savefig(dest / f'spatial_true_vs_estimated_seed{seed}.{ext}', dpi=220)
        plt.close(fig)
        a.file.close()
    df = pd.DataFrame(rows)
    df.to_csv(dest / 'metrics_by_seed.csv', index=False)
    df[['rmse', 'pearson_r', 'mean_absolute_error']].agg(['mean', 'std']).to_csv(dest / 'metrics_mean_sd.csv')
    print(df.to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
