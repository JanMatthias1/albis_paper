"""QC/HVG log-expression JSD across ten slices of the Figure 5A tissues."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import plot_distributions as base
OUT=base.FIGURE_5B_DIR/'jsd_overview_qc_hvg'

def sample_slice(sid):
    base.SLICE_ID=sid
    OUT.mkdir(exist_ok=True)
    arrays={};meta=[]
    for modality in ['cell','bin','spot']:
        data={}
        specs=[s for s in base.sources_for(modality) if modality=='cell' or s['label']!='splatter_sccube']
        for spec in specs:data[spec['label']]=base.load_and_qc(spec)
        # Same selection and transform order as the current QC/HVG composites.
        data=base.hvg_match_panels(data,cell_scale=modality=='cell')
        for spec in specs:
            a=data[spec['label']]
            if spec['label']!='splatter_sccube':
                sc.pp.normalize_total(a,target_sum=base.TARGET_SUM);sc.pp.log1p(a)
            values=a.X.data if hasattr(a.X,'tocsr') else np.asarray(a.X).ravel()
            values=values[values>0]
            assert values.size and np.isfinite(values).all()
            key=modality+'__'+spec['label']
            arrays[key]=base.cd.sample_values(values,base.SAMPLE_SIZE,np.random.default_rng(sid))
            meta.append(dict(modality=modality,method=spec['display_label'],n_obs=a.n_obs,n_genes=a.n_vars,source=str(spec['path'].resolve())))
        del data
    np.savez_compressed(OUT/f'slice_{sid}.npz',**arrays)
    (OUT/f'slice_{sid}.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('Finished slice',sid,flush=True)

def aggregate():
    for sid in range(10):
        if not (OUT/f'slice_{sid}.npz').exists():raise RuntimeError(f'Missing slice {sid}')
    upper=0
    for sid in range(10):
        with np.load(OUT/f'slice_{sid}.npz') as z:upper=max(upper,max(float(z[k].max()) for k in z.files))
    bins=np.linspace(0,upper,61);rows=[]
    for sid in range(10):
        with np.load(OUT/f'slice_{sid}.npz') as z:
            for modality in ['cell','bin','spot']:
                for method in (['spider','sccube'] if modality=='cell' else ['spider']):
                    value=base.cd.compute_jsd(z[modality+'__our_method'],z[modality+'__splatter_'+method],bins)
                    rows.append(dict(slice_id=sid,slice_number=sid+1,modality=modality,comparison='ALBIS vs '+('SPIDER' if method=='spider' else 'scCube'),jsd=value))
    df=pd.DataFrame(rows);assert len(df)==40 and df.jsd.between(0,1).all()
    df.to_csv(OUT/'jsd_by_slice.csv',index=False)
    fig,ax=plt.subplots(figsize=(8,5))
    for x,modality in enumerate(['cell','bin','spot']):
        for comparison,offset,color in [('ALBIS vs SPIDER',-.15,'#d95f02'),('ALBIS vs scCube',.15,'#7570b3')]:
            sub=df[(df.modality==modality)&(df.comparison==comparison)].sort_values('slice_id')
            if sub.empty:continue
            pos=x+offset
            bp=ax.boxplot([sub.jsd.to_numpy()],positions=[pos],widths=.24,patch_artist=True,showfliers=False,manage_ticks=False)
            bp['boxes'][0].set(facecolor=color,alpha=.3)
            ax.scatter(pos+np.linspace(-.065,.065,len(sub)),sub.jsd,color=color,s=25,zorder=3)
    ax.set(xticks=[0,1,2],xticklabels=['Cell','Bin','Spot'],ylabel='Jensen–Shannon divergence',ylim=(0,1.05),title='Log-normalized expression · QC / 556 genes')
    ax.legend(handles=[Patch(facecolor=c,alpha=.6,label=l) for l,c in [('ALBIS vs SPIDER','#d95f02'),('ALBIS vs scCube (cell only)','#7570b3')]],frameon=False,loc='best')
    fig.text(.5,.015,'10 slices per method and modality; points are slices',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.04,1,1))
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'jsd_boxplot.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    protocol={'n_slices':10,'scores':40,'replication':'Slices from one realization per method/modality; ALBIS uses separately calibrated Figure 2 tissues. Not independent simulation replicates.','qc':'Remove explicitly empty/unassigned observations.','genes':'556 per method; per-slice variance ranking as in current QC/HVG composite plots.','transform':'ALBIS/SPIDER: select genes using normalized-log variance, then normalize selected counts to 10000 and log1p, matching current QC/HVG composites. scCube: native reconstruction, ranked by native variance, no further transform.','histogram_edges':bins.tolist(),'sampling':'Nonzero values only; up to 2 million per method/slice/modality; seed=slice_id.','JSD':'Base-2 divergence (squared scipy Jensen-Shannon distance), shared histogram edges across all groups.','comparability':'scCube bin/spot excluded: sums of reconstructed cell log-expression.','note':'Shared bins and positive-entry sampling can change scores slightly relative to individually binned Slice 5 figures.'}
    (OUT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    print('Validated 40 scores; saved',OUT/'jsd_boxplot.png')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--slice-id',type=int,choices=range(10));p.add_argument('--aggregate',action='store_true');a=p.parse_args()
    if a.aggregate:aggregate()
    elif a.slice_id is not None:sample_slice(a.slice_id)
    else:p.error('Specify --slice-id or --aggregate')
