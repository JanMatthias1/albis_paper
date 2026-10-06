"""Combined spatial QC overview, verified against saved post-QC observations."""
import argparse
import json
from pathlib import Path
import sys
import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from anndata._io.specs import read_elem
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from manuscript_style import apply_style
apply_style()
ROOT=Path(__file__).resolve().parents[2]/'data/real_data_qc'
DATASETS=[('non_diseased_lung','Xenium · Non-diseased lung'),
 ('lung_cancer','Xenium · Lung cancer'),
 ('breast_cancer_visium_hd_16um','Visium HD · Breast cancer (16 µm)'),
 ('human_pancreas_visium_hd_16um','Visium HD · Pancreas (16 µm)'),
 ('lymph_node_visium','Visium · Lymph node'),('tonsil_visium','Visium · Tonsil')]
COLORS={'Retained':'#B9C0C7','Below median − 4 MAD':'#4C8FD5','Other Xenium exclusions':'#E65F5C','SpotSweeper exclusion':'#4C8FD5'}


def load_panel(name):
    folder=ROOT/name;summary=json.loads((folder/'qc_summary.json').read_text())
    files=list(folder.glob('*_qc.h5ad'));assert len(files)==1
    with h5py.File(files[0]) as f:
        kept=read_elem(f['obs']).index.astype(str);n_genes=len(read_elem(f['var']))
    raw=Path(summary['input_dir'])
    if 'nmads' in summary:
        import xenium_qc as xenium
        a=xenium.load_xenium_slice(raw)
        a,thresholds=xenium.apply_qc_filters(a,nmads=summary['nmads'])
        assert set(a.obs_names[~a.obs.low_qc])==set(kept),name
        obs=a.obs;xy=obs[['x_centroid','y_centroid']].to_numpy()
        categories=np.full(len(obs),'Retained',dtype=object)
        categories[obs.low_qc.to_numpy()]='Other Xenium exclusions'
        categories[(obs.low_total_count|obs.low_detected_features).to_numpy()]='Below median − 4 MAD'
        rule=f"Lower-tail {summary['nmads']:g}-MAD + controls"
    else:
        # The pre-QC universe is the filtered 10x barcode set restricted to tissue.
        matrix=list(raw.glob('*filtered_feature_bc_matrix.h5'));assert len(matrix)==1
        with h5py.File(matrix[0]) as f:
            barcodes=pd.Index(f['matrix/barcodes'].asstr()[:])
        spatial=raw/'spatial'
        if (spatial/'tissue_positions.parquet').exists():
            obs=pd.read_parquet(spatial/'tissue_positions.parquet').set_index('barcode')
        elif (spatial/'tissue_positions.csv').exists():
            obs=pd.read_csv(spatial/'tissue_positions.csv',index_col='barcode')
        else:
            obs=pd.read_csv(spatial/'tissue_positions_list.csv',header=None,index_col=0,
                names=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres'])
        obs=obs.loc[barcodes];obs=obs[obs.in_tissue==1]
        xy=obs[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy()
        assert set(kept).issubset(set(obs.index)),name
        categories=np.where(obs.index.isin(kept),'Retained','SpotSweeper exclusion')
        rule=f"SpotSweeper · k={summary['n_neighbors']}, cutoff={summary['cutoff']:g}"
    before=next(v for k,v in summary.items() if k.endswith('_before_qc'))
    assert len(xy)==before and int((categories=='Retained').sum())==len(kept),name
    return xy,categories,rule,dict(dataset=name,n_genes=n_genes,before_qc=before,after_qc=len(kept),
                                  excluded=before-len(kept),retained_percent=100*len(kept)/before)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--include-8um',action='store_true',default=True);args=ap.parse_args()
    datasets=DATASETS.copy()
    if args.include_8um:
        datasets += [('breast_cancer_visium_hd','Visium HD · Breast cancer (8 µm)'),('human_pancreas_visium_hd','Visium HD · Pancreas (8 µm)')]
    rows,cols=(3,3) if args.include_8um else (2,4)
    fig,axes=plt.subplots(rows,cols,figsize=(20,11) if rows==2 else (17,17))
    records=[]
    for i,((name,title),ax) in enumerate(zip(datasets,axes.flat)):
        xy,category,rule,record=load_panel(name);records.append(record)
        print(record,flush=True)
        # Small, equal-size markers and mixed drawing order avoid painting
        # every excluded observation over the dense retained population.
        rng = np.random.default_rng(2025)
        order = rng.permutation(len(xy))
        colors = np.array([COLORS[label] for label in category])
        size = 0.15 if '_hd' in name else (0.35 if 'visium' not in name else 9)
        ax.scatter(xy[order,0], xy[order,1], s=size, c=colors[order],
                   linewidths=0, rasterized=True, alpha=0.9)
        ax.set_aspect('equal');ax.invert_yaxis();ax.set_axis_off()
        ax.set_title(title,fontsize=17,pad=12)
        ax.text(.5,-.035,f"{record['n_genes']:,} genes\n{record['after_qc']:,} / {record['before_qc']:,} retained ({record['retained_percent']:.1f}%)\n{rule}",
                transform=ax.transAxes,ha='center',va='top',fontsize=14)
    handles=[Line2D([],[],marker='o',linestyle='',color=color,label=label,markersize=12) for label,color in COLORS.items()]
    # The legend fills the grid's first unused panel.
    key=axes.flat[len(datasets)];key.set_axis_off()
    key.legend(handles=handles,loc='center',frameon=False,fontsize=14,title='Spatial QC status',title_fontsize=16)
    for ax in axes.flat[len(datasets)+1:]: ax.set_axis_off()
    fig.subplots_adjust(left=.025,right=.985,top=.94,bottom=.06,wspace=.18,hspace=.40)
    out=ROOT/'supplementary';out.mkdir(exist_ok=True)
    plt.rcParams['svg.fonttype']='none'
    for ext in ['png','pdf','svg']:
        fig.savefig(out/f'real_data_qc_overview.{ext}',dpi=500,facecolor='white')
    plt.close(fig)
    pd.DataFrame(records).to_csv(out/'real_data_qc_overview_counts.csv',index=False)
    (out/'real_data_qc_overview_caption.txt').write_text(
        'Spatial quality-control overview of real-data sections. Gray observations were retained. '
        'Blue Xenium cells fall below median minus four unscaled median absolute deviations in log1p total counts or detected genes. '
        'Coral cells fail other Xenium criteria (zero count density or positive control-probe/codeword counts); blue takes precedence when criteria overlap. '
        'Blue Visium/Visium HD observations were excluded by SpotSweeper (36 neighbors; cutoff 3), '
        'using low total counts, low detected genes, or high mitochondrial percentage. '
        'Membership was verified against saved post-QC files. Maps have independent spatial scales and show every observation in the pre-QC tissue universe. Observations are drawn in a reproducibly shuffled order with equal marker sizes within each panel to reduce occlusion bias. '
        'Panel annotations report the number of genes and the retained counts and percentages.\n')
    print('Saved',out/'real_data_qc_overview.png',flush=True)


if __name__=='__main__': main()
