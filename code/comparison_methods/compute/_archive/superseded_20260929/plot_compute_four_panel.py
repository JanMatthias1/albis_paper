"""Four panels from fresh, consistently sampled single-seed benchmark processes."""
import argparse
import fcntl
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/compute-four-panel-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--task-only', action='store_true', help='Refresh available results after a task completes.')
args = parser.parse_args()
root = args.root.resolve()
out = root/'figures/four_panel_with_setup'
out.mkdir(parents=True, exist_ok=True)
with (root/'.report.lock').open('w') as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    records = [json.loads(p.read_text()) for p in sorted((root/'raw').glob('*/measurement.json'))]
    if not records:
        raise SystemExit('No completed measurements yet')
    data = pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,dict)} for r in records]).sort_values(['method','n_cells'])
    data.to_csv(root/'measurements.csv',index=False)
    methods = {'albis':('ALBIS','#8E63C7','o'), 'sccube':('scCube','#43B7A5','s'), 'spider':('SPIDER','#F2A65A','^')}
    panels = [
        ('simulation_seconds','Simulation time','VAE training and Splatter excluded','Elapsed time (s)','simulation_time'),
        ('total_seconds','Time including setup','Splatter + input + VAE training, where required','Elapsed time (s)','time_including_setup'),
        ('generation_peak_rss_gib','Memory during generation','Required reference and trained VAE remain resident','Peak RAM (GiB)','memory_during_generation'),
        ('including_setup_peak_rss_gib','Memory including setup','Peak across reference, setup and generation','Peak RAM (GiB)','memory_including_setup')]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'svg.fonttype':'none','pdf.fonttype':42})
    fig, axes = plt.subplots(2,2,figsize=(11,9))
    fig.subplots_adjust(left=.10,right=.97,bottom=.17,top=.84,wspace=.33,hspace=.85)
    sizes = json.loads((root/'protocol.json').read_text())['sizes']
    ticks = [10000,100000,600000,1000000,2000000,5000000]
    labels = ['10k','100k','600k','1M','2M','5M']
    for ax,(key,title,subtitle,ylabel,name) in zip(axes.flat,panels):
        ax.text(.5,1.22,title,transform=ax.transAxes,ha='center',fontsize=12,fontweight='bold')
        ax.text(.5,1.10,subtitle,transform=ax.transAxes,ha='center',fontsize=8,color='#555555')
        for method,(label,color,marker) in methods.items():
            # Reindex to show genuine gaps while jobs are pending.
            d=data.loc[data.method.eq(method)].set_index('n_cells').reindex(sizes)
            ax.plot(sizes,d[key],label=label,color=color,marker=marker,lw=1.7,markersize=5)
        ax.set_xscale('log')
        ax.set_xticks(ticks,labels,rotation=35)
        ax.set_xlim(8500,6000000)
        ax.set_xlabel('Requested cells (log scale)')
        ax.set_ylabel(ylabel)
        if key.endswith('seconds'):
            ax.set_yscale('log')
        else:
            ax.set_ylim(bottom=0)
        ax.grid(axis='y',color='#DDDDDD',lw=.6)
        ax.legend(frameon=False,fontsize=8)
    fig.text(.5,.96,'Computational cost of native tissue simulation',ha='center',fontsize=15,fontweight='bold')
    fig.text(.5,.92,f'{len(data)}/24 completed · CPU, 1 computational thread · 556 genes · seed 2025',ha='center',color='#555555')
    notes=('Fresh process for each method and size; scCube trains a fresh in-memory VAE. No saved models.\n'
           'Both memory panels use process RSS sampled every 10 ms; brief peaks may be missed.\n'
           'Setup-inclusive memory is the maximum across sequential Splatter and method phases, never their sum.\n'
           'Imports, reference export and validation are excluded. Generation retains required setup allocations.\n'
           'Native output representations differ. Single seed; no error bars. Missing points are pending runs.')
    fig.text(.07,.025,notes,fontsize=8,linespacing=1.5)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    bounds=[]
    for ax in axes.flat:
        headings=ax.texts[:2]
        for t in headings: t.set_visible(False)
        body=ax.get_tightbbox(renderer)
        center=(body.x0+body.x1)/2
        for t in headings:
            t.set_x(ax.transAxes.inverted().transform((center,body.y0))[0])
            t.set_visible(True)
        full=Bbox.union([body]+[t.get_window_extent(renderer) for t in headings])
        width=max(center-full.x0,full.x1-center)
        bounds.append(Bbox.from_extents(center-width,full.y0,center+width,full.y1).transformed(fig.dpi_scale_trans.inverted()).expanded(1.045,1.055))
    for ext in ['png','pdf','svg']:
        fig.savefig(out/f'compute_four_panel.{ext}',bbox_inches='tight',facecolor='white')
        for text in fig.texts:
            text.set_visible(False)
        for panel,bound in zip(panels,bounds):
            fig.savefig(out/f'{panel[-1]}.{ext}',bbox_inches=bound,facecolor='white')
        for text in fig.texts:
            text.set_visible(True)
    plt.close(fig)
    (out/'README.md').write_text(f'# Four-panel compute benchmark\n\n{len(data)}/24 completed.\n\n'+notes+'\n')
    (root/'RESULTS.md').write_text(f'# Fresh-process memory benchmark\n\n{len(data)}/24 completed.\n\n'+data[['method','n_cells','simulation_seconds','total_seconds','generation_peak_rss_gib','including_setup_peak_rss_gib']].to_csv(index=False)+'\n'+notes+'\n')
    print(f'{len(data)}/24 completed; {out}',flush=True)
