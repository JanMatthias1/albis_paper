"""Prepare the fixed 15-run, three-simulation-seed cell-type experiment."""
from pathlib import Path
import json,csv
ROOT=Path(__file__).resolve().parents[5]
PAPER=ROOT/'sim_paper'
OUT=PAPER/'data/figure_5/clustering/cell_type/experiment_600k'

def main():
    if OUT.exists():raise FileExistsError(OUT)
    OUT.mkdir()
    baseline=json.loads((PAPER/'data/figure_5/figure_5A_600k/settings.json').read_text())
    albis=json.loads((PAPER/'data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x/config.json').read_text())
    pools=[];runs=[]
    for seed in [2025,101,202]:
        for genes in [556,2000]:
            pool=OUT/f'pools/genes{genes}_seed{seed}';pool.mkdir(parents=True)
            contract={'n_cells':10000,'n_genes':genes,'n_cell_types':8,'cell_type_proportions':[.125]*8}
            (pool/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
            pools.append({'task':len(pools),'seed':seed,'genes':genes,'directory':str(pool)})
            for method in ['sccube','spider']:
                run=OUT/f'runs/{method}_genes{genes}_seed{seed}';run.mkdir(parents=True)
                cfg={k:baseline[k] for k in ['n_cells','n_slices','extent_um','bin_width_um','spot_spacing_um','spot_radius_um','sccube_grid_size','sccube_delta','sccube_lamda','sccube_epochs','sccube_cells_per_bin','sccube_cells_per_spot','spider_self_probability','spider_max_iterations']}
                cfg.update(seed=seed,n_genes=genes,proportions=[.125]*8,cell_only=True)
                (run/'settings.json').write_text(json.dumps(cfg,indent=2)+'\n')
                (run/'splatter').symlink_to(pool/'splatter',target_is_directory=True)
                runs.append({'task':len(runs),'method':method,'genes':genes,'seed':seed,'directory':str(run)})
        run=OUT/f'runs/albis_genes556_seed{seed}';run.mkdir(parents=True)
        cfg=dict(albis,batch_sigma=1.5,seed=seed,output_modalities=['cell'])
        (run/'settings.json').write_text(json.dumps(cfg,indent=2)+'\n')
        runs.append({'task':len(runs),'method':'albis','genes':556,'seed':seed,'directory':str(run)})
    for name,rows in [('pools',pools),('runs',runs)]:
        (OUT/f'{name}.json').write_text(json.dumps(rows,indent=2)+'\n')
    protocol={'n_cells':600000,'n_slices':10,'seeds':[2025,101,202],'pool_cells':10000,
        'conditions':['albis_556','sccube_556','sccube_2000','spider_556','spider_2000'],
        'albis':'Figure 5A uncropped 16um reference configuration; only batch_sigma=1.5, simulation seed and cell-only output selection changed.',
        'spider':'Native simulate_10X_3d; diagonal target .7, off-diagonal .3/7; no gyrus.',
        'analysis':'All cells; positive total counts; normalize10000/log1p/scale10/PCA30/Harmony(slice_id)/15NN/Leiden known-k=8; clustering seed 0.',
        'replication':'Simulation seeds; each scCube VAE trained independently. Shared Splatter pool between methods per gene-count/seed.',
        'limitations':['ALBIS has injected batch_sigma=1.5; competitors have no equivalent imposed batch effect.',
            'Gene count changes regenerate expression rather than subset the same matrix.',
            'SPIDER expression profiles are sampled from a 10k reference pool.']}
    (OUT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    print(OUT)
if __name__=='__main__':main()
