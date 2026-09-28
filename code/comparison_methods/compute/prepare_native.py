"""Write explicit native-call inputs and scheduler task lists; executes no methods."""
import argparse
import hashlib
import json
import math
from pathlib import Path

parser=argparse.ArgumentParser()
parser.add_argument('--root',type=Path,required=True)
parser.add_argument('--smoke',action='store_true')
args=parser.parse_args()
root=args.root.resolve();root.mkdir(parents=True,exist_ok=False)
for name in ['settings','raw','references','logs','code']:
    (root/name).mkdir()
code=Path(__file__).resolve().parent
project=code.parents[3]
base=json.loads((project/'sim_paper/data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x/config.json').read_text())
files=['albis_native.py','sccube_native.py','spider_native.py','reference_native.R',
       'albis_native.sbatch','sccube_native.sbatch','spider_native.sbatch','reference_native.sbatch',
       'prepare_native.py','submit_native.sh','report_native.py','report_native.sbatch']
hashes={}
for name in files:
    source=code/name;(root/'code'/name).write_bytes(source.read_bytes())
    hashes[name]=hashlib.sha256(source.read_bytes()).hexdigest()
protocol=dict(version='direct_native_v2',smoke=args.smoke,n_genes=556,
    reference_cells=128 if args.smoke else 10000,source_hashes=hashes,
    albis_source_sha256=hashlib.sha256((project/'albis/albis/simulation_sphere.py').read_bytes()).hexdigest(),
    output_policy='Native outputs unchanged. Missing native outputs unavailable. No composition reconstruction.',
    timing='Direct API call elapsed time; combined training/generation retained. Slurm adds whole-job accounting.',
    memory='Whole-process ru_maxrss; not isolated stage memory and not a process-tree sum.',
    coordinate_policy='ALBIS/SPIDER physical dimensions grow with N; scCube uses unchanged native coordinates/grid 8. No claim of matched physical density for scCube.',
    expression_policy='SPIDER native expression view is preserved; no forced copy/materialization.',
    scientific_scope='Costs of native representations, not identical outputs. No tissue export in timing.',
    threads=1,hardware_constraint='sapphirerapids')
(root/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
rows=['method\tn_cells\tseed\tsettings\toutput\treference']
sizes=[32] if args.smoke else [10000,50000,100000,200000,500000,600000,1000000]
seeds=[2025] if args.smoke else [2025,101,202]
for seed in seeds:
    for n in sizes:
        factor=(n/600000)**(1/3)
        albis=dict(base,n_cells=n,seed=seed,sphere_R_um=2050*factor,core_fuzz_width_um=102.5*factor,
                   batch_sigma=0,max_deg=0,max_shift=0,output_modalities=['bin','spot'],
                   capture_window_um=False,xenium_capture_window_um=False)
        settings=dict(protocol_version='direct_native_v2',n_cells=n,seed=seed,reference_cells=protocol['reference_cells'],
            sccube_epochs=1 if args.smoke else 200,spider_iterations=5 if args.smoke else 80000,
            target_cells_per_type={f'Group{i}':n//8 for i in range(1,9)},proportions=[.125]*8,
            spider_transition=[[.7 if i==j else .3/7 for j in range(8)] for i in range(8)],
            extent_um=(4*math.pi/3*2050**3*n/600000)**(1/3),albis=albis)
        path=root/'settings'/f'n{n}_seed{seed}.json';path.write_text(json.dumps(settings,indent=2)+'\n')
        for method in ['albis','sccube','spider']:
            rows.append('\t'.join(map(str,[method,n,seed,path,root/'raw'/f'{method}_n{n}_seed{seed}',root/'references'/f'seed{seed}'])))
(root/'tasks.tsv').write_text('\n'.join(rows)+'\n')
print(root)
