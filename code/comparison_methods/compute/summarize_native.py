#!/usr/bin/env python3
"""Collect native measurements, preserve failures, plot results, and review the pilot."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
from native_worker import save

METHODS = {'albis':('ALBIS · mRNA instances','#2c7fb8'),
           'sccube':('scCube · cell expression','#7570b3'),
           'spider':('SPIDER · cell expression','#d95f02')}


def accounting(root,submissions):
    ids = [r['job_id'] for r in submissions]
    if not ids: return {}
    result = subprocess.run(['sacct','-j',','.join(ids),'--parsable2','--noheader',
        '--format=JobID,State,ExitCode,ElapsedRaw,MaxRSS,NodeList,ReqMem,AllocCPUS'],capture_output=True,text=True)
    if result.returncode:
        (root/'accounting_error.log').write_text(result.stderr)
        return {}
    (root/'slurm_accounting.psv').write_text(result.stdout)
    parsed={}
    for line in result.stdout.splitlines():
        fields=line.split('|')
        if len(fields)>=8 and fields[0] in ids:
            parsed[fields[0]]=dict(zip(['job_id','state','exit_code','elapsed_seconds','max_rss','node','requested_memory','cpus'],fields[:8]))
    return parsed


def collect(root,fetch_accounting=True):
    cfg=json.loads((root/'protocol.json').read_text())
    tasks=json.loads((root/'tasks.json').read_text())
    submissions=json.loads((root/'submissions.json').read_text()) if (root/'submissions.json').exists() else []
    jobs={r['key']:r for r in submissions if r['kind']=='task'}
    scheduler=accounting(root,submissions) if fetch_accounting else {}
    rows=[]; stages=[]
    for task in tasks:
        directory=root/'raw'/task['key']
        path=directory/'worker/measurement.json'
        measured=json.loads(path.read_text()) if path.exists() else {}
        proc_path=directory/'process_measurement.json'
        process=json.loads(proc_path.read_text()) if proc_path.exists() else {}
        status=measured.get('status','submitted' if task['key'] in jobs else 'not_submitted')
        job=jobs.get(task['key'],{})
        slurm=scheduler.get(job.get('job_id'),{})
        state=slurm.get('state','').split()[0].rstrip('+') if slurm.get('state') else ''
        if status!='ok' and state in ['OUT_OF_MEMORY','TIMEOUT','CANCELLED','NODE_FAIL','FAILED','PREEMPTED']:
            status=state.lower()
        if status=='ok' and (process.get('returncode') not in [0,None] or state in ['OUT_OF_MEMORY','TIMEOUT','CANCELLED','NODE_FAIL','FAILED']):
            status='worker_endpoint_ok_job_failed'
        ok=status=='ok' and measured.get('endpoint_validated') and process.get('returncode')==0
        if status=='ok' and not ok: status='awaiting_process_completion'
        generation=[s for s in measured.get('stages',[]) if s['category']=='generation']
        pipeline=[s for s in measured.get('stages',[]) if s['category'] in ['input','generation','startup']]
        ref=measured.get('reference_measurement',{})
        output=measured.get('output',{})
        row=dict(key=task['key'],method=task['method'],n_cells=task['n_cells'],seed=task['seed'],
            protocol_id=cfg['protocol_id'],status=status,failure_reason=measured.get('failure_reason',''),
            generation_seconds=sum(s.get('wall_seconds',0) for s in generation) if ok else None,
            native_pipeline_seconds=sum(s.get('wall_seconds',0) for s in pipeline) if ok else None,
            reference_seconds=ref.get('wall_seconds',0) if ok else None,
            from_scratch_seconds=(sum(s.get('wall_seconds',0) for s in pipeline)+ref.get('wall_seconds',0)) if ok else None,
            generation_cpu_seconds=sum(s.get('cpu_seconds',0) for s in generation) if ok else None,
            generation_peak_gib=max([s['peak_tree_rss_bytes'] for s in generation] or [0])/2**30 if ok else None,
            from_scratch_peak_gib=max([s['peak_tree_rss_bytes'] for s in pipeline]+[ref.get('peak_tree_rss_bytes',0)])/2**30 if ok else None,
            process_wall_seconds=process.get('wall_seconds'), process_peak_gib=process.get('peak_tree_rss_bytes',0)/2**30 if process else None,
            n_genes=output.get('n_genes'),n_molecules=output.get('n_molecules'),actual_cells=output.get('n_cells'),
            expression_nnz=output.get('expression_nnz'),host=measured.get('host'),cpu_model=measured.get('cpu_model'),
            reference_checksum=ref.get('checksums',{}).get('counts.mtx'),job_id=job.get('job_id'),
            scheduler_state=state,source=str(path))
        rows.append(row)
        for stage in measured.get('stages',[]):
            stages.append(dict(key=task['key'],method=task['method'],n_cells=task['n_cells'],seed=task['seed'],**stage))
    return cfg,rows,stages


def csv_write(path,rows):
    if not rows: return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader();w.writerows(rows)


def review_pilot(root,cfg,rows):
    pilot=[r for r in rows if r['seed']==2025 and r['n_cells'] in [10000,100000]]
    result=dict(status='incomplete',resources={},warnings=[],basis='pilot observations with conservative scaling; not established complexity')
    if cfg['smoke_test']: return
    if len(pilot)!=6 or any(r['status']!='ok' for r in pilot):
        result['reason']='All six validated pilot runs must finish successfully.'
    else:
        result['status']='passed'
        cpus={r['cpu_model'] for r in pilot}
        if len(cpus)>1:
            result['warnings'].append('Pilot used multiple CPU models; runtime comparisons are confounded. Pin a hardware class before the full sweep.')
            result['status']='needs_review'
        for method in METHODS:
            small=next(r for r in pilot if r['method']==method and r['n_cells']==10000)
            large=next(r for r in pilot if r['method']==method and r['n_cells']==100000)
            exponent=max(1.,math.log(max(large['process_wall_seconds']/small['process_wall_seconds'],1),10))
            result['resources'][method]={}
            for n in cfg['sizes']:
                factor=max(1,n/100000)
                mem=4*math.ceil(max(16,large['process_peak_gib']*factor*2+4)/4)
                seconds=max(3600,large['process_wall_seconds']*factor**exponent*3)
                hours=math.ceil(seconds/3600)
                result['resources'][method][str(n)]=dict(memory=f'{mem}G',time_limit=f'{hours}:00:00',
                    extrapolation_exponent=exponent)
                if mem>512 or hours>72:
                    result['status']='needs_review'
                    result['warnings'].append(f'{method} at {n}: projected request {mem}G / {hours}h needs a separate resource decision.')
    save(root/'pilot_review.json',result)
    lines=['# Pilot review',f"Status: {result['status']}",'',result.get('reason','All six pilot endpoints validated.'),'']
    lines+=result['warnings']
    lines+=['','Resources in pilot_review.json are provisional requests, not measured requirements at larger sizes.']
    (root/'pilot_review.md').write_text('\n'.join(lines)+'\n')


def plot(root,rows,stages):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    done=[r for r in rows if r['status']=='ok']
    if not done: return
    plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none','font.size':10})
    fig,axes=plt.subplots(1,3,figsize=(14,4.4))
    metrics=[('generation_seconds','Native generation (seconds)'),('from_scratch_seconds','From scratch (seconds)'),
             ('from_scratch_peak_gib','Peak memory through native endpoint (GiB)')]
    for ax,(metric,label) in zip(axes,metrics):
        for method,(name,color) in METHODS.items():
            subset=[r for r in done if r['method']==method]
            sizes=sorted({r['n_cells'] for r in subset})
            means=[]; deviations=[]
            for n in sizes:
                values=[r[metric] for r in subset if r['n_cells']==n]
                ax.scatter([n]*len(values),values,color=color,s=14,alpha=.5)
                means.append(np.mean(values));deviations.append(np.std(values,ddof=1) if len(values)>1 else 0)
            ax.errorbar(sizes,means,yerr=deviations,color=color,marker='o',markersize=4,label=name,capsize=3)
        ax.set(xscale='log',yscale='log',xlabel='Requested tissue cells',ylabel=label)
        ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False)
    axes[0].legend(fontsize=8)
    fig.text(.5,.015,'Native outputs differ in resolution. scCube training included; reference preparation included in from-scratch cost.\nPoints show successful seeds; error bars are SD when replicated. Failures and success counts are in the CSV/report.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.1,1,1))
    for ext in ['png','pdf','svg']: fig.savefig(root/'figures'/f'native_scaling.{ext}',dpi=250)
    plt.close(fig)
    # Stage breakdown uses a common observed size, prioritizing the paper anchor.
    common=set.intersection(*[{r['n_cells'] for r in done if r['method']==method} for method in METHODS])
    if not common: return
    n=600000 if 600000 in common else max(common)
    selected=[s for s in stages if s['n_cells']==n and s.get('status')=='ok' and s['category'] in ['generation','input','startup']]
    names=list(dict.fromkeys(s['name'] for s in selected))
    fig,ax=plt.subplots(figsize=(9,4.5))
    bottoms=np.zeros(3)
    for name in names:
        values=[np.mean([s['wall_seconds'] for s in selected if s['method']==method and s['name']==name] or [0]) for method in METHODS]
        ax.bar(list(METHODS),values,bottom=bottoms,label=name.replace('_',' ')); bottoms+=values
    refs=[np.mean([r['reference_seconds'] for r in done if r['method']==method and r['n_cells']==n] or [0]) for method in METHODS]
    ax.bar(list(METHODS),refs,bottom=bottoms,label='reference generation')
    ax.set(ylabel='Mean wall time (seconds)',title=f'Native tissue preparation: {n:,} cells')
    ax.legend(bbox_to_anchor=(1.02,1),loc='upper left',fontsize=8);fig.tight_layout()
    for ext in ['png','pdf','svg']: fig.savefig(root/'figures'/f'native_stages_n{n}.{ext}',dpi=250)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--no-accounting',action='store_true')
    a=p.parse_args(); root=a.root.resolve()
    (root/'figures').mkdir(exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR',str(root/'cache/plots'))
    cfg,rows,stages=collect(root,not a.no_accounting)
    csv_write(root/'measurements.csv',rows);csv_write(root/'stages.csv',stages)
    review_pilot(root,cfg,rows);plot(root,rows,stages)
    lines=['# Native pre-aggregation benchmark results','',f"Protocol: {cfg['protocol_id']}",'',
        '| Method | Requested cells | Successful / submitted |','|---|---:|---:|']
    for method in METHODS:
        for n in cfg['sizes']:
            subset=[r for r in rows if r['method']==method and r['n_cells']==n and r['status']!='not_submitted']
            if subset: lines.append(f"| {method} | {n} | {sum(r['status']=='ok' for r in subset)} / {len(subset)} |")
    (root/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'successful':sum(r['status']=='ok' for r in rows),'total_planned':len(rows),'root':str(root)}))

if __name__=='__main__': main()
