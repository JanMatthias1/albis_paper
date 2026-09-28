#!/usr/bin/env python3
"""Native in-memory tissue endpoint. Run in the selected method environment."""
import argparse
from contextlib import contextmanager
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import random
import resource
import sys
import threading
import time
import traceback

ROOT = Path(__file__).resolve().parents[4]


def save(path, value):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


class Measurements:
    """Absolute RSS sum across this process tree, sampled at 50 ms and boundaries."""
    def __init__(self, out, report):
        self.pid = os.getpid()
        self.page_size = os.sysconf("SC_PAGE_SIZE")
        self.out, self.report = out, report
        self.active = None
        self.lock = threading.Lock()
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self.monitor, daemon=True)
        self.thread.start()

    def sample(self):
        pending, seen, rss = [self.pid], set(), 0
        while pending:
            pid = pending.pop()
            if pid in seen:
                continue
            seen.add(pid)
            try:
                proc = Path('/proc')/str(pid)
                rss += int((proc/'statm').read_text().split()[1]) * self.page_size
                for path in (proc/'task').glob('*/children'):
                    pending.extend(map(int,path.read_text().split()))
            except (OSError, ValueError, IndexError):
                pass
        with self.lock:
            if self.active is not None:
                self.active['peak_tree_rss_bytes'] = max(self.active['peak_tree_rss_bytes'], rss)
                self.active['rss_samples'] += 1

    def monitor(self):
        while not self.stop.wait(.05):
            self.sample()

    @staticmethod
    def cpu():
        a = resource.getrusage(resource.RUSAGE_SELF)
        b = resource.getrusage(resource.RUSAGE_CHILDREN)
        return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime

    @contextmanager
    def stage(self, name, category):
        stage = dict(name=name, category=category, status='running', peak_tree_rss_bytes=0, rss_samples=0)
        self.report['stages'].append(stage)
        save(self.out/'measurement.json', self.report)
        with self.lock:
            self.active = stage
        self.sample()
        cpu, started = self.cpu(), time.monotonic()
        try:
            yield
        except BaseException:
            stage['status'] = 'failed'
            raise
        else:
            stage['status'] = 'ok'
        finally:
            stage['wall_seconds'] = time.monotonic() - started
            stage['cpu_seconds'] = self.cpu() - cpu
            self.sample()
            with self.lock:
                self.active = None
            save(self.out/'measurement.json', self.report)

    def close(self):
        self.stop.set()
        self.thread.join()


def forbidden(*args, **kwargs):
    raise AssertionError('Pre-aggregation boundary crossed: capture/sectioning/containment called')


def albis(cfg, m, report):
    with m.stage('imports', 'startup'):
        import numpy as np
        sys.path.insert(0, str(ROOT/'albis'))
        import albis.simulation_sphere as sim
    report['source_hashes'][str(Path(sim.__file__).resolve())] = digest(sim.__file__)
    for name in ['section_3d_molecule_sphere', 'assign_points_to_cells_by_containment',
                 'aggregate_molecules_to_cell_level', 'aggregate_molecules_to_grid_bins_2d_slices_window',
                 'aggregate_molecules_to_spots_2d_slices_window']:
        setattr(sim, name, forbidden)
    params = dict(cfg['albis'])
    factor = (cfg['n_cells']/600000)**(1/3)
    params.update(n_cells=cfg['n_cells'], seed=cfg['seed'], sphere_R_um=2050*factor,
                  core_fuzz_width_um=102.5*factor, batch_sigma=0, max_deg=0, max_shift=0,
                  capture_window_um=False, xenium_capture_window_um=False,
                  output_modalities=['bin', 'spot'])
    report['resolved_method_settings'] = params
    with m.stage('native_molecular_tissue', 'generation'):
        result = sim.simulate_3d_molecule_sphere_base(**params)
    with m.stage('validate', 'validation'):
        cells = result['adata_cell_true']
        mol = result['molecules']
        xyz = mol['full_xyz']
        assert cells.shape == (cfg['n_cells'], 556)
        assert xyz.shape == (len(mol['full_gene']), 3)
        assert len(xyz) == int(cells.X.sum()) == result['meta']['n_molecules_generated_total']
        assert len(mol['assigned_gene']) == 0
        assert result['meta']['n_molecules_assigned_total'] is None
        for start in range(0, len(xyz), 100000):
            end = start + 100000
            assert np.isfinite(xyz[start:end]).all()
            assert ((mol['full_gene'][start:end] >= 0) & (mol['full_gene'][start:end] < 556)).all()
        assert np.isfinite(cells.obsm['spatial']).all()
        assert np.unique(cells.obsm['spatial'][:,2]).size > min(10,cfg['n_cells']-1)
        assert (cells.X.data >= 0).all() and np.isfinite(cells.X.data).all()
        report['output'] = dict(n_cells=cells.n_obs, n_genes=cells.n_vars, n_molecules=len(xyz),
            expression_nnz=int(cells.X.nnz), expression_dtype=str(cells.X.dtype),
            native_representation='explicit_mrna_instances', coordinates_shape=list(xyz.shape),
            radius_um=params['sphere_R_um'], full_stream_array_bytes=sum(a.nbytes for k,a in mol.items() if k.startswith('full_')),
            aggregation_called=False)


def competitor(method, cfg, reference_dir, m, report):
    with m.stage('imports', 'startup'):
        import numpy as np
        import pandas as pd
        import anndata as ad
        from scipy.io import mmread
        if method == 'spider':
            import spider
            import spider.sim_expr as expr_module
            assert importlib.metadata.version('st-spider') == '1.2.0'
        else:
            import torch
            from scCube.sccube import scCube
            sccube_module = importlib.import_module('scCube.sccube')
    with m.stage('reference_loading_preprocessing', 'input'):
        X = mmread(reference_dir/'counts.mtx').T.tocsr()
        genes = (reference_dir/'genes.tsv').read_text().splitlines()
        meta = pd.read_csv(reference_dir/'cells.tsv', sep='\t')
        types = [f'type{i}' for i in range(1,9)]
        assert X.shape == (cfg['reference_cells'], 556)
        assert len(genes) == 556 and len(set(genes)) == 556
        assert set(meta.cell_type_true) == set(types)
        if method == 'spider':
            reference = ad.AnnData(X, obs=meta.set_index('cell_id'), var=pd.DataFrame(index=genes))
        else:
            torch.manual_seed(cfg['seed'])
            torch.set_num_threads(1)
            torch.set_num_interop_threads(1)
            model = scCube()
            model.generate_spot_data_random = forbidden
            reference = model.pre_process(pd.DataFrame(X.T.toarray(), index=genes, columns=meta.cell_id),
                pd.DataFrame({'Cell':meta.cell_id, 'Cell_type':meta.cell_type_true}), is_normalized=False)
    extent = (4*np.pi/3 * 2050**3 * cfg['n_cells']/600000)**(1/3)
    if method == 'spider':
        import inspect
        for func in [spider.simulate_10X_3d, spider.get_sim_cell_level_expr]:
            source = inspect.getsourcefile(func)
            report['source_hashes'][source] = digest(source)
        spider.get_sim_spot_level_expr = forbidden
        expr_module.get_sim_spot_level_expr = forbidden
        prior = np.full(8, 1/8)
        trans = np.full((8,8), .3/7)
        np.fill_diagonal(trans, .7)
        report['resolved_method_settings'] = dict(extent_um=extent, prior=prior.tolist(), target_trans=trans.tolist(),
            smallsample_max_iter=cfg['spider_iterations'], bigsample_max_iter=10000, tol=.02, T=.01)
        with m.stage('native_3d_placement', 'generation'):
            labels, xyz = spider.simulate_10X_3d(cell_num=cfg['n_cells'], Num_celltype=8, prior=prior,
                target_trans=trans, image_width=extent, image_height=extent, image_depth=extent,
                smallsample_max_iter=cfg['spider_iterations'])
        with m.stage('native_expression_sampling_materialization', 'generation'):
            expression = spider.get_sim_cell_level_expr(celltype_assignment=labels, adata=reference,
                Num_celltype=8, Num_ct_sample=np.bincount(labels,minlength=8), match_list=types,
                ct_key='cell_type_true').copy()
            values = expression.X
            labels = np.asarray(types)[labels]
    else:
        report['source_hashes'][sccube_module.__file__] = digest(sccube_module.__file__)
        import inspect
        source = inspect.getsourcefile(sccube_module.generate_vae)
        report['source_hashes'][source] = digest(source)
        sizes = np.full(8, cfg['n_cells']//8, dtype=int)
        sizes[:cfg['n_cells'] % 8] += 1
        report['resolved_method_settings'] = dict(extent_um=extent, epochs=cfg['sccube_epochs'],
            spatial_dim=3, spatial_size=8, delta=2.0, lamda=.75, is_split=False, used_device='cpu',
            target_num=dict(zip(types,map(int,sizes))), save_model=False)
        with m.stage('native_training_and_expression', 'generation'):
            generated_meta, data = model.train_vae_and_generate_cell(reference,
                celltype_key='Cell_type',cell_key='Cell', target_num=dict(zip(types,map(int,sizes))),
                epoch_num=cfg['sccube_epochs'], used_device='cpu', save_model=False)
        with m.stage('native_3d_placement', 'generation'):
            data, generated_meta = model.generate_pattern_random(data, generated_meta,
                spatial_dim=3, spatial_size=8, delta=2., lamda=.75, is_split=False, set_seed=True, seed=cfg['seed'])
            generated_meta = generated_meta.loc[data.columns]
            xyz = generated_meta[['point_x','point_y','point_z']].to_numpy() * (extent/8)
            values = data.T.to_numpy()
            labels = generated_meta.Cell_type.to_numpy()
    with m.stage('validate', 'validation'):
        assert values.shape == (cfg['n_cells'],556)
        assert xyz.shape == (cfg['n_cells'],3) and np.isfinite(xyz).all()
        assert np.unique(xyz[:,2]).size > min(10,cfg['n_cells']-1)
        assert len(labels) == cfg['n_cells'] and set(labels).issubset(set(types))
        if method == 'sccube':
            assert data.columns.is_unique and generated_meta.index.equals(data.columns)
            assert list(data.index) == genes
        nnz = 0
        for start in range(0,cfg['n_cells'],10000):
            block = values[start:start+10000]
            v = block.data if hasattr(block,'tocsr') else block
            assert np.isfinite(v).all() and (v>=0).all()
            nnz += int(np.count_nonzero(v))
        if method == 'spider':
            assert np.array_equal(expression.obs.cell_type_true.to_numpy(), labels)
        report['output'] = dict(n_cells=cfg['n_cells'], n_genes=556,n_molecules=None,
            expression_nnz=nnz, expression_dtype=str(values.dtype),native_representation='cell_level_expression',
            coordinates_shape=list(xyz.shape), extent_um=extent,aggregation_called=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--settings',type=Path,required=True)
    p.add_argument('--method',choices=['albis','sccube','spider'],required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--reference',type=Path)
    a = p.parse_args()
    cfg = json.loads(a.settings.read_text())
    a.out.mkdir(parents=True,exist_ok=False)
    report = dict(status='running',method=a.method,n_cells=cfg['n_cells'],seed=cfg['seed'],
        protocol_id=cfg['protocol_id'],settings=cfg,stages=[],source_hashes={str(Path(__file__).resolve()):digest(__file__)},
        memory_measurement='sum of process-tree RSS, 50ms samples and stage boundaries; shared pages may count twice',
        host=platform.node(),cpu_model=next((x.split(':',1)[1].strip() for x in Path('/proc/cpuinfo').read_text().splitlines() if x.startswith('model name')),'unknown'),
        slurm={k:v for k,v in os.environ.items() if k in ['SLURM_JOB_ID','SLURM_CPUS_PER_TASK','SLURM_MEM_PER_NODE','SLURM_JOB_PARTITION']},
        thread_settings={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']})
    m = Measurements(a.out,report)
    try:
        for package, expected in cfg.get('environment_versions', {}).get(a.method, {}).items():
            assert importlib.metadata.version(package) == expected, f'Package version changed: {package}'
        with m.stage('seed_initialization', 'startup'):
            import numpy as np
            random.seed(cfg['seed']); np.random.seed(cfg['seed'])
        if a.method == 'albis':
            albis(cfg,m,report)
        else:
            ref = json.loads((a.reference/'reference_measurement.json').read_text())
            assert ref['status'] == 'ok' and ref['seed'] == cfg['seed']
            assert ref['n_genes'] == 556 and ref['n_cells'] == cfg['reference_cells']
            # Validate inputs before timing generation; record checksum verification separately.
            with m.stage('reference_integrity', 'startup'):
                for name, expected in ref['checksums'].items():
                    assert digest(a.reference/name) == expected, name
            report['reference_measurement'] = ref
            competitor(a.method,cfg,a.reference,m,report)
        report['status'] = 'ok'
        report['endpoint_validated'] = True
    except BaseException as e:
        report['status'] = 'failed'; report['failure_reason'] = f'{type(e).__name__}: {e}'
        traceback.print_exc()
    finally:
        m.close()
        report['versions'] = {}
        for name in ['numpy','scipy','anndata','torch','scCube','st-spider','psutil']:
            try: report['versions'][name] = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError: pass
        save(a.out/'measurement.json',report)
    if report['status'] != 'ok':
        sys.exit(1)
    print(json.dumps(report['output']),flush=True)


if __name__ == '__main__':
    main()
