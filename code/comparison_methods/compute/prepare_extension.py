"""Prepare single-seed 200k/400k/600k native inputs and frozen source snapshots."""
import hashlib
import json
from pathlib import Path
import shutil
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--validated-split', type=Path, required=True)
args = parser.parse_args()
root, source = args.root.resolve(), args.source.resolve()
validation = json.loads((args.validated_split / 'validation.json').read_text())
assert validation['status'] == 'passed'
root.mkdir(parents=True, exist_ok=False)
for folder in ['settings', 'raw', 'logs', 'code']:
    (root / folder).mkdir()
code = Path(__file__).resolve().parent
files = ['albis_native.py', 'spider_native.py', 'sccube_split.py',
         'native_extension.sbatch', 'sccube_extension.sbatch', 'prepare_extension.py',
         'report_native_extension.py', 'report_native_extension.sbatch',
         'plot_compute_three_panel.py']
for name in files:
    shutil.copy2(code / name, root / 'code' / name)
template = json.loads((source / 'settings/n600000_seed2025.json').read_text())
rows = ['method\tn_cells\tseed\tsettings\toutput\treference']
for n in [200000, 400000, 600000]:
    s = json.loads(json.dumps(template))
    factor = (n / 600000) ** (1 / 3)
    s.update(n_cells=n, target_cells_per_type={f'Group{i}': n // 8 for i in range(1, 9)},
             extent_um=template['extent_um'] * factor)
    s['albis'].update(n_cells=n, sphere_R_um=template['albis']['sphere_R_um'] * factor,
                      core_fuzz_width_um=template['albis']['core_fuzz_width_um'] * factor)
    (root / 'settings' / f'n{n}_seed2025.json').write_text(json.dumps(s, indent=2) + '\n')
    for method in ['albis', 'spider']:
        rows.append('\t'.join(map(str, [method, n, 2025, root / 'settings' / f'n{n}_seed2025.json',
                    root / 'raw' / f'{method}_n{n}_seed2025', source / 'references/seed2025'])))
(root / 'tasks.tsv').write_text('\n'.join(rows) + '\n')
protocol = dict(seed=2025, sizes=[200000, 400000, 600000], n_genes=556,
                threads=1, hardware_constraint='sapphirerapids',
                source_pilot=str(source), validated_split=str(args.validated_split.resolve()),
                prior_validation=validation, sccube_epochs=200,
                sccube_policy='Train once; retain VAE in memory across new sizes; no saved models.',
                reference_policy='Reuse original seed2025 reference; original Splatter cost reported separately.',
                native_output_policy='Unchanged native outputs; SPIDER view retained; no aggregation/export.',
                memory_policy='scCube cumulative process RSS; ALBIS/SPIDER independent process RSS.',
                source_hashes={name: hashlib.sha256((root / 'code' / name).read_bytes()).hexdigest()
                               for name in files})
(root / 'protocol.json').write_text(json.dumps(protocol, indent=2) + '\n')
print(root)
