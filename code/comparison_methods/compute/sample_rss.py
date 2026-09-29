"""Read a worker's /proc RSS every 100 ms; no interaction with model objects."""
import argparse
import csv
import os
from pathlib import Path
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pid', type=int, required=True)
parser.add_argument('--phase-file', type=Path, required=True)
parser.add_argument('--done-file', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
page = os.sysconf('SC_PAGE_SIZE')
with args.output.open('w', buffering=1) as f:
    writer = csv.writer(f)
    writer.writerow(['monotonic_seconds', 'phase', 'rss_bytes'])
    while not args.done_file.exists():
        try:
            rss = int(Path(f'/proc/{args.pid}/statm').read_text().split()[1]) * page
            phase = args.phase_file.read_text().strip() if args.phase_file.exists() else 'imports_and_input'
            writer.writerow([time.monotonic(), phase or 'phase_transition', rss])
        except (FileNotFoundError, ProcessLookupError):
            break
        time.sleep(0.1)
