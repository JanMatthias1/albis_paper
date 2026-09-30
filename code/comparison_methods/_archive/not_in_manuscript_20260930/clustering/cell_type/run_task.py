"""Execute one manifest task with the method's established environment."""
import argparse,json,subprocess,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5]
CODE=Path(__file__).resolve().parent
OUT=ROOT/'sim_paper/data/figure_5/clustering/cell_type/experiment_600k'
def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['pools','runs']);p.add_argument('task',type=int);a=p.parse_args()
    row=json.loads((OUT/f'{a.stage}.json').read_text())[a.task];run=Path(row['directory'])
    env=os.environ.copy()
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:env[key]='4'
    env['PYTHONHASHSEED']='0'
    for key,sub in [('MPLCONFIGDIR','matplotlib'),('NUMBA_CACHE_DIR','numba'),('XDG_CACHE_HOME','xdg')]:
        path=run/'cache'/sub;path.mkdir(parents=True,exist_ok=True);env[key]=str(path)
    def call(args,log):
        with (run/log).open('w') as f:subprocess.run(list(map(str,args)),env=env,stdout=f,stderr=subprocess.STDOUT,check=True,cwd=ROOT)
    if a.stage=='pools':
        call([ROOT/'comparison_methods/env/splatter/bin/Rscript',ROOT/'comparison_methods/code/workflows/splatter_expression.R','--contract',run/'contract.json','--seed',row['seed'],'--out-dir',run/'splatter'],'generation.log')
    else:
        method=row['method'];python=ROOT/'sim_paper/env/albis-tutorial/bin/python'
        if method=='albis':call([python,'-u',CODE/'generate_albis.py','--run',run],'generation.log')
        else:call([ROOT/f'comparison_methods/env/{method}/bin/python','-u',CODE.parent.parent/'overview/generate.py','--method',method,'--settings',run/'settings.json','--out',run],'generation.log')
        call([python,'-u',CODE/'cluster.py','--run',run,'--method',method,'--genes',row['genes']],'clustering.log')
    print('Completed',row,flush=True)
if __name__=='__main__':main()
