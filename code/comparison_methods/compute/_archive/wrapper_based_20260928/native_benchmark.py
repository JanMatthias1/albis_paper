#!/usr/bin/env python3
"""Prepare, submit, or execute the versioned native pre-aggregation benchmark."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from native_worker import save, digest

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[3]
BENCH = ROOT/'comparison_methods'
PYTHON = BENCH/'env/analysis/bin/python'
ENVIRONMENTS = {'albis':'our_method','sccube':'sccube','spider':'spider'}


def environment(out, seed):
    env = dict(os.environ)
    env.update({name:'1' for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS']})
    env.update(PYTHONHASHSEED=str(seed), PYTHONDONTWRITEBYTECODE='1',
        MPLCONFIGDIR=str(out/'cache/matplotlib'), NUMBA_CACHE_DIR=str(out/'cache/numba'), XDG_CACHE_HOME=str(out/'cache'))
    for p in ['cache/matplotlib','cache/numba']:
        (out/p).mkdir(parents=True,exist_ok=True)
    return env


def prepare(root, smoke=False):
    root.mkdir(parents=True,exist_ok=False)
    for d in ['settings','raw','references','logs','code']:
        (root/d).mkdir()
    base = json.loads((ROOT/'sim_paper/data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x/config.json').read_text())
    protocol = dict(protocol_id='native_preaggregation_v1'+('_smoke' if smoke else ''),n_genes=556,
        reference_cells=128 if smoke else 10000,n_cell_types=8, computational_threads=1, hardware_constraint='sapphirerapids',
        sccube_epochs=1 if smoke else 200, spider_iterations=5 if smoke else 80000,
        albis=base, seeds=[2025] if smoke else [2025,101,202],
        sizes=[32] if smoke else [10000,50000,100000,200000,500000,600000,1000000],
        density_cells_per_um3=600000/(4*3.141592653589793/3*2050**3),
        output_policy='native in-memory endpoint; compact summaries only; no tissue export',
        primary_cost='generation including scCube training; input preparation and startup recorded separately',
        reference_policy='fixed synthetic reference per seed; charged once to each competitor from-scratch total',
        smoke_test=smoke, implementation_hashes={})
    for name in ['native_worker.py','native_benchmark.py','native_job.sh','summarize_native.py','advance_native.py','finish_native.sh']:
        source = CODE/name
        (root/'code'/name).write_bytes(source.read_bytes())
        protocol['implementation_hashes'][name] = digest(source)
    protocol['environment_versions'] = {}
    version_code = '''import importlib.metadata as m,json
result={}
for name in ['numpy','scipy','anndata','torch','scCube','st-spider']:
 try: result[name]=m.version(name)
 except m.PackageNotFoundError: pass
print(json.dumps(result))'''
    for method, envname in ENVIRONMENTS.items():
        protocol['environment_versions'][method] = json.loads(subprocess.check_output(
            [BENCH/'env'/envname/'bin/python','-c',version_code],text=True))
    # Code snapshots are audit copies. Execute canonical code only after hash checks.
    protocol['albis_source_sha256'] = digest(ROOT/'albis/albis/simulation_sphere.py')
    save(root/'protocol.json',protocol)
    tasks = []
    for seed in protocol['seeds']:
        for n in protocol['sizes']:
            settings = dict(protocol,n_cells=n,seed=seed)
            path = root/'settings'/f'n{n}_seed{seed}.json'
            save(path,settings)
            for method in ENVIRONMENTS:
                tasks.append(dict(method=method,n_cells=n,seed=seed,key=f'{method}_n{n}_seed{seed}',settings=str(path)))
    save(root/'tasks.json',tasks)
    print(root,flush=True)


def verify_code(root):
    protocol = json.loads((root/'protocol.json').read_text())
    for name, expected in protocol['implementation_hashes'].items():
        if digest(CODE/name) != expected:
            raise RuntimeError(f'Implementation changed after preparation: {name}; use a new protocol directory')
    if digest(ROOT/'albis/albis/simulation_sphere.py') != protocol['albis_source_sha256']:
        raise RuntimeError('ALBIS source changed after preparation')
    return protocol


def timed_command(cmd,out,env):
    import psutil
    started = time.monotonic()
    peak = 0
    with (out/'stdout.log').open('w') as stdout, (out/'stderr.log').open('w') as stderr:
        proc = subprocess.Popen(['/usr/bin/time','-v','-o',str(out/'time_v.txt'),*map(str,cmd)],
            stdout=stdout,stderr=stderr,env=env,cwd=ROOT)
        parent = psutil.Process(proc.pid)
        while proc.poll() is None:
            rss = 0
            try:
                processes = [parent] + parent.children(recursive=True)
            except psutil.Error:
                processes = []
            for p in processes:
                try: rss += p.memory_info().rss
                except psutil.Error: pass
            peak = max(peak,rss)
            time.sleep(.05)
        returncode = proc.wait()
    measured = dict(returncode=returncode,wall_seconds=time.monotonic()-started,peak_tree_rss_bytes=peak)
    # External high-water RSS is a secondary process measurement, not a tree sum.
    for line in (out/'time_v.txt').read_text().splitlines():
        if 'Maximum resident set size (kbytes):' in line:
            measured['gnu_time_max_rss_bytes'] = int(line.rsplit(':',1)[1])*1024
        for label,key in [('User time (seconds):','user_seconds'),('System time (seconds):','sys_seconds')]:
            if label in line: measured[key] = float(line.rsplit(':',1)[1])
    return measured


def reference(root,seed):
    cfg = verify_code(root)
    out = root/'references'/f'seed{seed}'
    out.mkdir(exist_ok=False)
    report = dict(status='running',seed=seed,n_cells=cfg['reference_cells'],n_genes=556,protocol_id=cfg['protocol_id'])
    save(out/'reference_measurement.json',report)
    contract = dict(n_cells=cfg['reference_cells'],n_genes=556,n_cell_types=8,cell_type_proportions=[.125]*8)
    save(out/'contract.json',contract)
    source = BENCH/'code/workflows/splatter_expression.R'
    report['source_sha256'] = digest(source)
    result = timed_command([BENCH/'env/splatter/bin/Rscript',source,'--contract',out/'contract.json',
        '--seed',seed,'--out-dir',out],out,environment(out,seed))
    report.update(result)
    manifest_path = out/'manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    ok = result['returncode']==0 and manifest.get('status')=='ok'
    if ok:
        from scipy.io import mmread
        import pandas as pd
        matrix = mmread(out/'counts.mtx')
        meta = pd.read_csv(out/'cells.tsv',sep='\t')
        genes = (out/'genes.tsv').read_text().splitlines()
        ok = matrix.shape==(556,cfg['reference_cells']) and len(set(genes))==556 and set(meta.cell_type_true)=={f'type{i}' for i in range(1,9)}
        report['checksums'] = {name:digest(out/name) for name in ['counts.mtx','genes.tsv','cells.tsv']}
    report['manifest'] = manifest
    report['status'] = 'ok' if ok else 'failed'
    save(out/'reference_measurement.json',report)
    if not ok: raise RuntimeError('Reference generation or validation failed')


def run_task(root,key):
    verify_code(root)
    task = next(t for t in json.loads((root/'tasks.json').read_text()) if t['key']==key)
    out = root/'raw'/key
    out.mkdir(exist_ok=False)
    save(out/'task.json',dict(task,status='running',job_id=os.environ.get('SLURM_JOB_ID')))
    cmd = [BENCH/'env'/ENVIRONMENTS[task['method']]/'bin/python',CODE/'native_worker.py',
        '--settings',task['settings'],'--method',task['method'],'--out',out/'worker']
    if task['method'] != 'albis': cmd += ['--reference',root/'references'/f"seed{task['seed']}"]
    result = timed_command(cmd,out,environment(out,task['seed']))
    save(out/'process_measurement.json',result)
    path = out/'worker/measurement.json'
    measurement = json.loads(path.read_text()) if path.exists() else {}
    ok = result['returncode']==0 and measurement.get('status')=='ok'
    save(out/'task.json',dict(task,status='ok' if ok else 'failed',job_id=os.environ.get('SLURM_JOB_ID')))
    if not ok: raise RuntimeError(f'Task failed: {key}; see {out}')


def submit(root, phase, methods=None):
    cfg = verify_code(root)
    tasks = json.loads((root/'tasks.json').read_text())
    if phase=='pilot': tasks = [t for t in tasks if t['seed']==2025 and t['n_cells'] in [10000,100000]]
    elif phase=='smoke':
        if not cfg['smoke_test']: raise ValueError('Not a smoke protocol')
    elif phase=='full':
        gate = root/'pilot_review.json'
        review = json.loads(gate.read_text()) if gate.exists() else {}
        if review.get('status')!='passed': raise RuntimeError('A successful pilot review is required before full submission')
    if methods: tasks = [t for t in tasks if t['method'] in methods]
    records_path = root/'submissions.json'
    records = json.loads(records_path.read_text()) if records_path.exists() else []
    def sbatch(kind,key,seed,dependency=None,mem='16G',hours='02:00:00'):
        name = f'native_{key}'
        cmd = ['sbatch','--parsable','--job-name',name,'--kill-on-invalid-dep=yes','--partition','shared','--cpus-per-task','2',
            '--constraint',cfg['hardware_constraint'],'--mem',mem,'--time',hours,'--output',str(root/'logs'/f'{key}_%j.out')]
        if dependency: cmd += ['--dependency',f'afterok:{dependency}']
        cmd += [str(CODE/'native_job.sh'),str(root),kind,str(key if kind=='task' else seed)]
        job = subprocess.check_output(cmd,text=True).strip().split(';')[0]
        if not job.isdigit(): raise RuntimeError(f'Unexpected sbatch response: {job}')
        records.append(dict(kind=kind,key=key,seed=seed,job_id=job,phase=phase,memory=mem,time_limit=hours,dependency=dependency))
        save(records_path,records)
        print(json.dumps(records[-1]),flush=True)
        return job
    pools = {r['seed']:r['job_id'] for r in records if r['kind']=='reference'}
    for seed in sorted({t['seed'] for t in tasks}):
        if seed not in pools: pools[seed] = sbatch('reference',f'reference_seed{seed}',seed)
    submitted = {r['key'] for r in records if r['kind']=='task'}
    # At most three tissue workers concurrently: serialize sizes within each method.
    previous = {method:None for method in ENVIRONMENTS}
    for task in tasks:
        if task['key'] in submitted: continue
        dependencies = []
        refpath = root/'references'/f"seed{task['seed']}"/'reference_measurement.json'
        reference_ready = refpath.exists() and json.loads(refpath.read_text()).get('status')=='ok'
        if task['method']!='albis' and not reference_ready:
            dependencies.append(pools[task['seed']])
        if previous[task['method']]: dependencies.append(previous[task['method']])
        if phase=='full':
            resources = json.loads((root/'pilot_review.json').read_text())['resources'][task['method']][str(task['n_cells'])]
            mem,hours = resources['memory'],resources['time_limit']
        else: mem,hours = ('16G','00:30:00') if phase=='smoke' else ('64G','1-00:00:00')
        previous[task['method']] = sbatch('task',task['key'],task['seed'],':'.join(dependencies) or None,mem,hours)
    if phase=='pilot' and not (root/'pilot_monitor.json').exists():
        ids=':'.join(r['job_id'] for r in records if r['phase']=='pilot' and r['kind']=='task')
        job=subprocess.check_output(['sbatch','--parsable','--job-name','native_pilot_review','--partition','shared',
            '--cpus-per-task','1','--mem','8G','--time','00:30:00','--dependency',f'afterany:{ids}',
            '--output',str(root/'logs/pilot_review_%j.out'),str(CODE/'finish_native.sh'),str(root),'pilot'],text=True).strip().split(';')[0]
        if not job.isdigit(): raise RuntimeError(job)
        save(root/'pilot_monitor.json',dict(job_id=job,action='review pilot; advance only if passed; schedule final reporting'))
        print(f'Pilot review and continuation job: {job}',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['prepare','reference','task','submit'])
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--smoke',action='store_true')
    p.add_argument('--seed',type=int,default=2025)
    p.add_argument('--key')
    p.add_argument('--method',nargs='+',choices=list(ENVIRONMENTS))
    p.add_argument('--phase',choices=['smoke','pilot','full'],default='pilot')
    a=p.parse_args(); root=a.root.resolve()
    if a.command=='prepare': prepare(root,a.smoke)
    elif a.command=='reference': reference(root,a.seed)
    elif a.command=='task': run_task(root,a.key)
    else: submit(root,a.phase,a.method)

if __name__=='__main__': main()
