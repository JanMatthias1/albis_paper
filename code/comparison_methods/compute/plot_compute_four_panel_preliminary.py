"""Preliminary four-panel preview from existing measurements; memory proxies explicitly labeled."""
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
out = root/'figures/four_panel_preliminary'
out.mkdir(parents=True, exist_ok=True)
with (root/'.report.lock').open('w') as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    data = pd.read_csv(root/'measurements.csv')
    data = data.loc[data.status.eq('ok')].sort_values(['method','n_cells'])
    data['total_seconds'] = data['first_use_seconds']
    data['generation_peak_rss_gib'] = data['generation_sampled_peak_rss_gib']
    # Whole-process peaks are explicit proxies, never presented as phase measurements.
    other = data.method.ne('sccube')
    data.loc[other, 'generation_peak_rss_gib'] = data.loc[other, 'process_peak_rss_gib']
    data['including_setup_peak_rss_gib'] = data['process_peak_rss_gib']
    data.to_csv(out/'plotted_measurements.csv',index=False)
    methods = {'albis':('ALBIS','#8E63C7','o'), 'sccube':('scCube','#43B7A5','s'), 'spider':('SPIDER','#F2A65A','^')}
    panels = [
        ('simulation_seconds','Simulation time','VAE training and Splatter excluded','Elapsed time (s)','simulation_time'),
        ('total_seconds','Time including setup','Splatter + input + VAE training, where required','Elapsed time (s)','time_including_setup'),
        ('generation_peak_rss_gib','Memory during generation','Hollow markers: whole-process proxy','Peak RAM (GiB)','memory_during_generation'),
        ('including_setup_peak_rss_gib','Memory including setup','Partial: method-process peak; Splatter unmeasured','Peak RAM (GiB)','memory_including_setup')]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'svg.fonttype':'none','pdf.fonttype':42})
    fig, axes = plt.subplots(2,2,figsize=(11,9))
    fig.subplots_adjust(left=.10,right=.97,bottom=.17,top=.84,wspace=.33,hspace=.85)
    sizes = sorted(data.n_cells.unique())
    ticks = [10000,100000,600000,1000000,2000000,5000000]
    labels = ['10k','100k','600k','1M','2M','5M']
    for ax,(key,title,subtitle,ylabel,name) in zip(axes.flat,panels):
        ax.text(.5,1.22,title,transform=ax.transAxes,ha='center',fontsize=12,fontweight='bold')
        ax.text(.5,1.10,subtitle,transform=ax.transAxes,ha='center',fontsize=8,color='#555555')
        for method,(label,color,marker) in methods.items():
            # Reindex to show genuine gaps while jobs are pending.
            d=data.loc[data.method.eq(method)].set_index('n_cells').reindex(sizes)
            memory = key.endswith('gib')
            proxy = memory and (key == 'including_setup_peak_rss_gib' or method != 'sccube')
            if key == 'generation_peak_rss_gib':
                label += ' (100 ms samples)' if method == 'sccube' else ' (process peak proxy)'
            elif key == 'including_setup_peak_rss_gib' and method == 'sccube':
                label += ' (includes batch history)'
            ax.plot(sizes,d[key],label=label,color=color,marker=marker,lw=1.7,markersize=5,
                    markerfacecolor='white' if proxy else color,
                    linestyle='--' if proxy else '-')
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
    fig.text(.5,.92,'PRELIMINARY · Existing single-seed results · CPU, 1 thread · 556 genes',ha='center',color='#555555')
    notes=('Timing panels use all 24 existing measurements. Setup time is charged at each size; VAE trained once per batch.\n'
           'Generation memory: scCube measured only at 1M–5M; smaller sizes are unavailable. ALBIS/SPIDER use process-peak proxies.\n'
           'Setup memory is partial: Splatter RSS was not measured. scCube includes training and earlier work in each batch.\n'
           'Memory panels are a provisional comparison, not matching phase-isolated measurements. No new runs or seeds.\n'
           'Native output representations differ. Hollow markers / dashed lines flag proxies or incomplete setup coverage.')
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
    (out/'README.md').write_text(f'# Preliminary four-panel compute benchmark\n\n{len(data)}/24 completed.\n\n'+notes+'\n')
    (out/'RESULTS.md').write_text(f'# Preliminary memory comparison\n\n{len(data)}/24 completed.\n\n'+data[['method','n_cells','simulation_seconds','total_seconds','generation_peak_rss_gib','including_setup_peak_rss_gib']].to_csv(index=False)+'\n'+notes+'\n')
    print(f'{len(data)}/24 completed; {out}',flush=True)
