"""Paths and JSON/task helpers for the isolated Figure 3 experiment."""
import csv
import importlib.util
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[4]
CODE = Path(__file__).resolve().parent
BASE = ROOT / 'albis_paper/data/figure_3/no_harmony_domain'
SHARED = ROOT / 'albis_paper/code/clustering'
SEEDS = [2025, 101, 202]
MODALITIES = ['cell', 'bin16um', 'spot']
K_GEOM = dict(cell=60, bin16um=100, spot=8)
LAMBDAS = [0., .1, .2, .3, .5, .8, 1.]

def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)

def load_shared(name, filename):
    import sys
    sys.path.insert(0, str(SHARED))
    sys.path.insert(0, str(SHARED.parent))
    spec = importlib.util.spec_from_file_location(name, SHARED / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def task_for(root, phase, index):
    rows = json.loads((root / f'{phase}_tasks.json').read_text())
    task = rows[index]
    assert task['task'] == index
    if phase == 'final' and task['pipeline'] == 'banksy':
        selected = json.loads((root / 'selected_parameters.json').read_text())
        task['lambda'] = selected['modalities'][task['modality']]['lambda']
    return task
