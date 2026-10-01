"""Collect all36 final results, flag cluster mismatches, and compare prior Harmony runs."""
import argparse
import json
from pathlib import Path
import pandas as pd
from common import *

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=BASE);root=p.parse_args().root
    rows=[];status=[]
    for task in json.loads((root/'final_tasks.json').read_text()):
        resolved=task_for(root,'final',task['task']);path=Path(task['directory'])/'metrics.json'
        state=dict(task=task['task'],directory=task['directory'],status='missing')
        if path.exists():
            m=json.loads(path.read_text());assert not m['harmony_applied']
            assert m['task']['input']==resolved['input'] and m['task']['lambda']==resolved['lambda']
            # not_converged: Leiden missed exactly 6 domains within the resolution search's 15 trials;
            # kept in metrics_by_seed.csv but excluded from the means and figures
            state['status']='complete' if m['achieved_clusters']==6 else 'not_converged'
            rows.append(dict(mix=task['mix'],modality=task['modality'],seed=task['seed'],pipeline=task['pipeline'],**{'lambda':resolved['lambda']},k_geom=task['k_geom'],n_obs=m['n_obs'],domain_ari=m['domain_ari'],slice_ari=m['slice_ari'],achieved_clusters=m['achieved_clusters'],status=state['status'],path=str(path)))
        status.append(state)
    pd.DataFrame(status).to_csv(root/'summary/status.csv',index=False)
    df=pd.DataFrame(rows);df.to_csv(root/'summary/metrics_by_seed.csv',index=False)
    if len(rows)!=36:raise RuntimeError('Final results incomplete; inspect summary/status.csv')
    df[df.status=='complete'].groupby(['mix','modality','pipeline']).agg(mean_ari=('domain_ari','mean'),sample_sd=('domain_ari','std'),n_seeds=('seed','size'),min_clusters=('achieved_clusters','min'),max_clusters=('achieved_clusters','max')).to_csv(root/'summary/metrics_mean_sd.csv')
    # Paired seed values at the old selected lambdas may differ from the new selection.
    import sys
    sys.path.insert(0,str(SHARED/'ari_recovery_summary'))
    import plot_banksy_vs_pca_recovery as previous
    old,old_lambdas=previous.load();comparison=[]
    for r in rows:
        key=(r['mix'],'banksy' if r['pipeline']=='banksy' else 'genes_bs0',r['modality'],'domain_true')
        value=old.get(key,{}).get(r['seed'])
        comparison.append(dict(r,previous_harmony_ari=value,previous_lambda=old_lambdas[(r['modality'],'domain_true')] if r['pipeline']=='banksy' else None,ari_difference=None if value is None else r['domain_ari']-value,interpretation='BANKSY parameters may also change after retuning; not an isolated Harmony effect.'))
    pd.DataFrame(comparison).to_csv(root/'summary/comparison_to_previous_harmony.csv',index=False)
    excluded=[r for r in rows if r['status']=='not_converged']
    for r in excluded:print(f"Not converged (excluded from means): {r['mix']}/{r['modality']}/seed{r['seed']}/{r['pipeline']}, {r['achieved_clusters']} clusters",flush=True)
    print(f'Validated 36 final runs; {36-len(excluded)} converged to 6 domains.',flush=True)

if __name__=='__main__':main()
