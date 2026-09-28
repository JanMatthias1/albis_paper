"""Read measurement files only; never loads, augments or transforms model outputs."""
import argparse
import csv
import json
from pathlib import Path

parser=argparse.ArgumentParser()
parser.add_argument('--root',type=Path,required=True)
args=parser.parse_args();root=args.root
with (root/'tasks.tsv').open() as f:
    tasks=list(csv.DictReader(f,delimiter='\t'))
submitted={}
for receipt in root.glob('submissions_*.tsv'):
    with receipt.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if row['kind']=='task': submitted[row['key']]=row['job_id']
rows=[]
for task in tasks:
    path=Path(task['output'])/'measurement.json'
    key=Path(task['output']).name
    record=json.loads(path.read_text()) if path.exists() else {}
    reference_path=Path(task['reference'])/'measurement.json'
    ref=json.loads(reference_path.read_text()) if reference_path.exists() else {}
    rows.append(dict(method=task['method'],n_cells=int(task['n_cells']),seed=int(task['seed']),
        job_id=submitted.get(key),status=record.get('status','submitted_or_failed' if key in submitted else 'not_submitted'),
        generation_seconds=record.get('generation_seconds'),input_seconds=record.get('input_seconds',0) if record else None,
        reference_generation_seconds=ref.get('generation_seconds') if task['method']!='albis' else 0,
        process_peak_rss_gib=record['process_peak_rss_bytes']/2**30 if record else None,
        n_genes=record.get('n_genes'),n_molecules=record.get('n_molecules'),
        endpoint=record.get('endpoint'),source=str(path)))
with (root/'measurements.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
(root/'RESULTS.md').write_text('# Direct-native compute run\n\n'
    +'Successful runs: '+str(sum(r['status']=='ok' for r in rows))+' / '+str(len(submitted)-sum(k.startswith('seed') for k in submitted))+' submitted tissue runs.\n\n'
    +'Missing measurements are not successes. Check scheduler/logs for pending or failed jobs. '
    +'Generation timings exclude input loading and imports; scCube training is included. '
    +'Reference generation is listed separately. Peak RSS is for the whole method process. '
    +'Native outputs remain unchanged and are not exported by this timing-only experiment.\n')
print(root/'measurements.csv')
