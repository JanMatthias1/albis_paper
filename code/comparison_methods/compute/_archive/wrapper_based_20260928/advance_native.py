#!/usr/bin/env python3
"""Continue only after successful pilot review; schedule final reporting."""
import argparse
import json
from pathlib import Path
import subprocess
from native_worker import save
from native_benchmark import CODE, submit


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();root=a.root.resolve()
    review=json.loads((root/'pilot_review.json').read_text())
    if review['status']!='passed':
        save(root/'continuation.json',dict(status='held_for_review',reason=review.get('reason'),warnings=review.get('warnings')))
        print('Full sweep held: inspect pilot_review.md',flush=True)
        return
    submit(root,'full')
    records=json.loads((root/'submissions.json').read_text())
    dependencies=':'.join(r['job_id'] for r in records if r['kind']=='task')
    result=subprocess.check_output(['sbatch','--parsable','--job-name','native_compute_report','--partition','shared',
        '--cpus-per-task','1','--mem','8G','--time','00:30:00','--dependency',f'afterany:{dependencies}',
        '--output',str(root/'logs/final_report_%j.out'),str(CODE/'finish_native.sh'),str(root),'final'],text=True)
    job=result.strip().split(';')[0]
    if not job.isdigit(): raise RuntimeError(result)
    save(root/'continuation.json',dict(status='full_submitted',final_report_job=job))
    print(f'Final report job: {job}',flush=True)

if __name__=='__main__': main()
