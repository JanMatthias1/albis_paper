#!/usr/bin/env python3
"""Short full-size native STAGATE memory check on the common Figure 4 tissue."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import runpy
import sys
import time
import traceback

CODE=Path(__file__).resolve().parent
ROOT=CODE.parents[3]
INPUT=ROOT/'sim_paper/data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x/data/cell/simulation_cell_z_qc.h5ad'


def save(path,report):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(report,indent=2)+'\n');tmp.replace(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--arm',choices=['3d','2d'],required=True)
    parser.add_argument('--epochs',type=int,default=5)
    a=parser.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    path=a.out/f'feasibility_{a.arm}.json'
    if path.exists(): raise FileExistsError(path)
    report=dict(status='running',input=str(INPUT),input_mtime_ns=INPUT.stat().st_mtime_ns,
        arm=a.arm,epochs=a.epochs,subsample=False,training_seed=0,graph_model='knn',k_2d=6,k_z=3,
        hidden_dims=[512,30],job_id=os.environ.get('SLURM_JOB_ID'),
        scope='Full-size graph, native training and native final inference; no mclust/UMAP or quality claim',
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),CODE/'3D_stagate.py']})
    save(path,report)
    started=time.monotonic()
    try:
        sys.path.insert(0,str(CODE))
        functions=runpy.run_path(str(CODE/'3D_stagate.py'))
        import scanpy as sc
        import torch
        import numpy as np
        if not torch.cuda.is_available(): raise RuntimeError('GPU required for this feasibility check')
        report['gpu_name']=torch.cuda.get_device_name(0)
        report['gpu_total_gib']=torch.cuda.get_device_properties(0).total_memory/2**30
        report['torch_version']=torch.__version__
        timer=time.monotonic()
        adata=sc.read_h5ad(INPUT)
        assert adata.shape==(600000,556)
        assert len(adata.obs.slice_id.unique())==10
        adata=functions['prepare_input'](adata)
        report['load_preprocess_seconds']=time.monotonic()-timer
        report['n_cells'],report['n_genes']=adata.shape
        timer=time.monotonic()
        functions['cal_spatial_net_3d_knn'](adata,k_2d=6,k_z=3)
        report['graph_seconds']=time.monotonic()-timer
        report['edges_3d']=len(adata.uns['Spatial_Net'])
        report['edges_2d']=len(adata.uns['Spatial_Net_2D'])
        graph=adata.uns['Spatial_Net' if a.arm=='3d' else 'Spatial_Net_2D']
        report['selected_edges']=len(graph)
        report['status']='training'
        save(path,report)
        torch.cuda.synchronize();torch.cuda.reset_peak_memory_stats()
        timer=time.monotonic()
        embedding=functions['run_stagate'](adata,graph,a.epochs,'STAGATE_check')
        torch.cuda.synchronize()
        report['native_training_and_inference_seconds']=time.monotonic()-timer
        report['gpu_peak_allocated_gib']=torch.cuda.max_memory_allocated()/2**30
        report['gpu_peak_reserved_gib']=torch.cuda.max_memory_reserved()/2**30
        assert embedding.shape==(600000,30) and np.isfinite(embedding).all()
        report['embedding_shape']=list(embedding.shape)
        report['status']='ok'
    except BaseException as error:
        report['status']='failed';report['error']=f'{type(error).__name__}: {error}'
        traceback.print_exc()
        try:
            import torch
            report['gpu_peak_allocated_gib']=torch.cuda.max_memory_allocated()/2**30
            report['gpu_peak_reserved_gib']=torch.cuda.max_memory_reserved()/2**30
        except Exception: pass
    finally:
        report['total_wall_seconds']=time.monotonic()-started
        report['host_peak_rss_gib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/2**30
        save(path,report)
    print(json.dumps(report,indent=2),flush=True)
    if report['status']!='ok': sys.exit(1)

if __name__=='__main__': main()
