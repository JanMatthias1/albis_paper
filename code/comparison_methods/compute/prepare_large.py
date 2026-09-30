"""Prepare 1M/2M/5M seed-2025 native benchmarking inputs and immutable worker snapshots."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--previous', type=Path, required=True)
args = parser.parse_args()
root, previous = args.root.resolve(), args.previous.resolve()
prior = json.loads((previous / 'protocol.json').read_text())
source = Path(prior['source_pilot'])
if not source.exists():
    source = previous.parent / source.name
validated_split = Path(prior['validated_split'])
if not validated_split.exists():
    validated_split = previous.parent / validated_split.name
root.mkdir(parents=True, exist_ok=False)
for folder in ['settings', 'raw', 'logs', 'code']:
    (root / folder).mkdir()
code = Path(__file__).resolve().parent
files = ['albis_native.py', 'spider_native.py', 'sccube_split.py', 'sample_rss.py',
         'native_extension.sbatch', 'sccube_large.sbatch', 'prepare_large.py',
         'report_large.py', 'report_large.sbatch', 'plot_compute_four_panel_preliminary.py']
for name in files:
    shutil.copy2(code / name, root / 'code' / name)
template = json.loads((source / 'settings/n600000_seed2025.json').read_text())
sizes = [1000000, 2000000, 5000000]
rows = ['method\tn_cells\tseed\tsettings\toutput\treference']
for n in sizes:
    s = json.loads(json.dumps(template))
    factor = (n / 600000) ** (1/3)
    s.update(n_cells=n, target_cells_per_type={f'Group{i}': n // 8 for i in range(1, 9)},
             extent_um=template['extent_um'] * factor)
    s['albis'].update(n_cells=n, sphere_R_um=template['albis']['sphere_R_um'] * factor,
                      core_fuzz_width_um=template['albis']['core_fuzz_width_um'] * factor)
    settings = root / 'settings' / f'n{n}_seed2025.json'
    settings.write_text(json.dumps(s, indent=2) + '\n')
    for method in ['albis', 'spider']:
        rows.append('\t'.join(map(str, [method, n, 2025, settings,
                    root / 'raw' / f'{method}_n{n}_seed2025', source / 'references/seed2025'])))
(root / 'tasks.tsv').write_text('\n'.join(rows) + '\n')
protocol = dict(seed=2025, sizes=sizes, n_genes=556, threads=1, hardware_constraint='sapphirerapids',
                previous_run=str(previous), source_pilot=str(source),
                validated_split=str(validated_split), sccube_epochs=200,
                scCube_memory='100-ms external sampling of worker RSS by phase; VAE resident. Also retain cumulative ru_maxrss. RSS includes process allocations retained from previous stages.',
                scCube_gate='Before 2M and 5M: prior high-water RSS × size ratio × 1.75 < 80% of 96 GiB allocation.',
                scheduling='ALBIS/SPIDER submitted one size at a time after reviewing prior results.',
                source_hashes={name: hashlib.sha256((root / 'code' / name).read_bytes()).hexdigest() for name in files},
                native_output_policy=prior['native_output_policy'], reference_policy=prior['reference_policy'])
(root / 'protocol.json').write_text(json.dumps(protocol, indent=2) + '\n')
print(root)
