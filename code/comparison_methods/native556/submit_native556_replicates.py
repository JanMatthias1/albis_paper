"""Submit prepared seed-101/202 pools, generation/RCTD, and clustering jobs."""
import json
from pathlib import Path
import subprocess

PAPER=Path(__file__).resolve().parents[3]
OVERVIEW=PAPER/'code/comparison_methods/native556'
CLUSTER=PAPER/'code/comparison_methods/clustering/cell_type'
manifest=PAPER/'data/comparison_methods/native556_replicate_jobs.json'
if manifest.exists():
    raise SystemExit(f'Submission record already exists: {manifest}; inspect before resubmitting')
records=[]
manifest.write_text('[]\n')
def submit(args, record):
    result=subprocess.run(['sbatch','--parsable',*map(str,args)],check=True,capture_output=True,text=True)
    job=result.stdout.strip().split(';')[0]
    if not job.isdigit(): raise RuntimeError(result.stdout)
    records.append(dict(record,job_id=job))
    manifest.write_text(json.dumps(records,indent=2)+'\n')
    print(records[-1],flush=True)
    return job
for seed in [101,202]:
    out=PAPER/f'data/comparison_methods/overview_native556_seed{seed}'
    assert json.loads((out/'settings.json').read_text())['seed']==seed
    pool=submit([f'--output={OVERVIEW}/logs/native556_seed{seed}_pool_%j.out',
                 OVERVIEW/'submit_native556_pool.sh',out,seed],dict(seed=seed,stage='pool'))
    for method in ['spider','sccube']:
        generation=submit([f'--dependency=afterok:{pool}',
            f'--export=ALL,NATIVE556_RUN_PREFIX=native556_seed{seed}',
            f'--output={OVERVIEW}/logs/native556_seed{seed}_{method}_%j.out',
            OVERVIEW/'submit_native556_method.sh',out,method,seed],
            dict(seed=seed,method=method,stage='generation_rctd_plots'))
        submit([f'--dependency=afterany:{generation}',
            f'--export=ALL,CELL_CLUSTER_SOURCE={out},CELL_CLUSTER_PREFIX=native556_seed{seed},CELL_CLUSTER_N_GENES=556',
            f'--output={CLUSTER}/logs/native556_seed{seed}_{method}_%j.out',
            CLUSTER/'run_pilot.sh',method],dict(seed=seed,method=method,stage='clustering'))
