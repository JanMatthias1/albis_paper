#!/usr/bin/env python
"""Cell-only Figure 5B comparison on the native scCube training scale."""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import plot_distributions as base


def main():
    specs=base.sources_for('cell')
    datasets=[]
    for spec in specs:
        a=base.load_and_qc(spec)
        if spec['label']=='splatter_sccube':
            transform='Native VAE reconstruction of log1p(normalize_total(10000)); unchanged'
        else:
            # Normalize on the full native panel BEFORE any gene selection.
            sc.pp.normalize_total(a,target_sum=base.TARGET_SUM)
            sc.pp.log1p(a)
            transform='normalize_total(10000) on full native panel, then log1p'
        mean,var=base.cd.gene_mean_var(a.X)
        assert np.isfinite(mean).all() and np.isfinite(var).all()
        datasets.append(dict(spec,mean=mean,var=var,genes=a.var_names.to_numpy(),n_obs=a.n_obs,transform=transform))
        del a
    target=min(len(d['genes']) for d in datasets)
    for matched in [False,True]:
        folder=base.FIGURE_5B_DIR/'cell'/('distribution_comparison_qc_hvg' if matched else 'distribution_comparison')
        folder.mkdir(parents=True,exist_ok=True)
        fig,ax=plt.subplots(figsize=(8,6.5))
        rows=[];meta=[]
        for d in datasets:
            idx=np.arange(len(d['genes']))
            if matched and len(idx)>target:idx=np.sort(np.argsort(d['var'],kind='stable')[-target:])
            mean,var=d['mean'][idx],d['var'][idx]
            # Linear axes: the data themselves are already on a log-expression scale.
            ax.scatter(mean,var,s=10,alpha=.45,linewidths=0,color=d['color'],label=d['display_label'],rasterized=True)
            rows.append(pd.DataFrame({'method':d['display_label'],'gene':d['genes'][idx],'mean_log_normalized_expression':mean,'variance_log_normalized_expression':var}))
            meta.append({'method':d['display_label'],'input':str(d['path'].resolve()),'n_cells':d['n_obs'],'native_genes':len(d['genes']),'plotted_genes':len(idx),'transform':d['transform']})
        ax.set(xlabel='Mean log-normalized expression',ylabel='Variance of log-normalized expression',title='Cell expression mean–variance · Slice 5',xlim=(0,None),ylim=(0,None))
        ax.legend(frameon=False);fig.tight_layout()
        for ext in ['png','pdf','svg']:fig.savefig(folder/f'mean_variance_log_normalized.{ext}',dpi=220,bbox_inches='tight')
        plt.close(fig)
        pd.concat(rows,ignore_index=True).to_csv(folder/'mean_variance_log_normalized_genes.csv',index=False)
        summary={'slice_id':base.SLICE_ID,'slice_number':base.SLICE_ID+1,'modality':'cell','matched_panel':matched,'sources':meta,'reference_curves':None,'axes':'linear; expression values are already log-transformed','selection':'Top variance on the comparison scale, after transformation; no re-normalization after selection' if matched else 'All native genes','caveats':['scCube is a VAE reconstruction on the training scale, not observed normalized counts.','Gene identities and native panel sizes differ; matching gene count does not match biology.','Matched scCube gene selection differs from the historical raw-panel plot: no second normalization/log transform.','Not applied to bins/spots because scCube aggregates already transformed cell values.']}
        (folder/'mean_variance_log_normalized_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print('Saved',folder/'mean_variance_log_normalized.png',flush=True)
if __name__=='__main__':main()
