"""Figure 5B: distributions on the log-normalized expression scale only."""
import json
from pathlib import Path
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import plot_distributions as base

def main():
    panels=[]
    for modality in ['cell','bin','spot']:
        data=[]
        for spec in base.sources_for(modality):
            if modality!='cell' and spec['label']=='splatter_sccube':continue
            a=base.load_and_qc(spec)
            reconstructed=spec['label']=='splatter_sccube'
            if not reconstructed:
                sc.pp.normalize_total(a,target_sum=base.TARGET_SUM);sc.pp.log1p(a)
            _,var=base.cd.gene_mean_var(a.X)
            data.append((spec,a,var,reconstructed))
        target=min(a.n_vars for _,a,_,_ in data)
        for matched in [False,True]:
            samples=[]
            for spec,a,var,reconstructed in data:
                idx=np.sort(np.argsort(var,kind='stable')[-target:]) if matched and a.n_vars>target else np.arange(a.n_vars)
                X=a.X[:,idx]
                values=X.data if hasattr(X,'tocsr') else np.asarray(X).ravel()
                values=values[values>0]
                sample=base.cd.sample_values(values,base.SAMPLE_SIZE,np.random.default_rng(0))
                samples.append(dict(spec,values=sample,n_obs=a.n_obs,n_genes=len(idx),reconstructed=reconstructed))
            panels.append((modality,matched,samples))
        del data
    # Identical histogram edges for every modality/panel, including zero as origin.
    edges=np.linspace(0,max(float(d['values'].max()) for _,_,ds in panels for d in ds),61)
    for modality,matched,ds in panels:
        folder=base.FIGURE_5B_DIR/modality/('distribution_comparison_qc_hvg' if matched else 'distribution_comparison')
        fig,ax=plt.subplots(figsize=(9,6))
        for d in ds:
            label=d['display_label']+(' (VAE reconstruction)' if d['reconstructed'] else '')
            ax.hist(d['values'],bins=edges,density=True,color=d['color'],alpha=.5,label=label)
        pairs=[(d['display_label'],base.cd.compute_jsd(ds[0]['values'],d['values'],edges)) for d in ds[1:]]
        base.annotate_jsd_pairs(ax,pairs,loc='upper right')
        ax.set(title=f'Log-normalized expression · {modality.capitalize()} · Slice 5',xlabel='Expression value (nonzero entries)',ylabel='Density')
        ax.legend(frameon=False,loc='upper left');fig.tight_layout()
        for ext in ['png','pdf','svg']:fig.savefig(folder/f'log_normalized_expression_compare.{ext}',dpi=200,bbox_inches='tight')
        plt.close(fig)
        meta={'modality':modality,'slice_id':4,'matched_panel':matched,'histogram_edges':edges.tolist(),'JSD_vs_ALBIS':dict(pairs),'transform':'ALBIS/SPIDER normalize_total(10000) on full native panel, then log1p. scCube cell reconstruction unchanged. Gene selection after transformation; no renormalization.','sampling':'Positive entries only; up to 2 million per method, random seed 0.','scCube_aggregate':'Excluded: sums of cell log-expression are not log-normalized aggregate counts.','sources':[{'method':d['display_label'],'path':str(d['path'].resolve()),'n_obs':d['n_obs'],'n_genes':d['n_genes']} for d in ds]}
        (folder/'log_normalized_expression_summary.json').write_text(json.dumps(meta,indent=2)+'\n')
        old=folder/'raw_norm_log_compare.png'
        if old.exists():
            archive=folder/'archive';archive.mkdir(exist_ok=True)
            old.rename(archive/'raw_norm_log_compare_before_log_only.png')
        print('Saved',folder/'log_normalized_expression_compare.png',flush=True)
if __name__=='__main__':main()
