#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
from pathlib import Path
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    from tqdm import tqdm
except Exception:
    tqdm = None

parser = argparse.ArgumentParser(description="Sum TimeSlotBusy (s) per sample from .err/.log pairs")
parser.add_argument("dir", nargs="?", default=".", help="directory with .err and .log files (default: .)")
parser.add_argument("-j", "--jobs", type=int, default=28, help="number of worker threads for I/O (default: 8)")
parser.add_argument("--tail-bytes", type=int, default=65536, help="how many bytes from file end to search for TimeSlotBusy (default: 65536)")
parser.add_argument("--no-progress", action="store_true", help="disable progress bar")
args = parser.parse_args()

d = Path(args.dir)
if not d.is_dir():
    print("Not a directory:", d, file=sys.stderr)
    sys.exit(2)

err_files = sorted(d.glob("*.err"))

re_sample = re.compile(r'\bsample\b\s+(\S+)', re.IGNORECASE)
re_time = re.compile(r'TimeSlotBusy\s*\(s\)\s*:\s*([0-9]+)')
# pattern to find finished time in .err files, e.g. "Plots finished in 558.19s" or "finished in 558.19 seconds"
re_err_time = re.compile(r'finished in\s*([0-9]+(?:\.[0-9]+)?)\s*(?:s|seconds?)', re.IGNORECASE)


def extract_seconds_from_err_text(text: str) -> float:
    """Extract the last "finished in ... s/seconds" value from given text.
    Returns seconds as float, or 0.0 if not found.
    """
    try:
        matches = re_err_time.findall(text)
        if matches:
            return float(matches[-1])
    except Exception:
        return 0.0
    return 0.0


def process_err(err_path: Path):
    """Parse .err to get sample name and extract finished seconds from the .err text.
    Returns (sample, time_seconds) or None on failure.
    """
    try:
        text = err_path.read_text(errors="ignore")
    except Exception:
        return None
    m = re_sample.search(text)
    if not m:
        return None
    sample = m.group(1)
    # extract finished time from the .err text itself
    time_val = extract_seconds_from_err_text(text)
    # return seconds as float
    return sample, time_val


results = []
workers = max(1, args.jobs)
use_progress = (tqdm is not None) and (not args.no_progress)

if workers == 1:
    iterator = err_files
    if use_progress:
        iterator = tqdm(err_files, desc="Processing .err files")
    for err in iterator:
        r = process_err(err)
        if r:
            results.append(r)
else:
    futures = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for err in err_files:
            futures.append(ex.submit(process_err, err))
        if use_progress:
            it = tqdm(as_completed(futures), total=len(futures), desc="Processing .err files")
            for fut in it:
                r = fut.result()
                if r:
                    results.append(r)
        else:
            for fut in as_completed(futures):
                r = fut.result()
                if r:
                    results.append(r)

# aggregate values per sample (collect seconds)
sample_values = defaultdict(list)
for s, t in results:
    sample_values[s].append(t)

# group samples with specified prefixes
group_prefixes = [
    'EGamma_', 'Muon_', 'MuonEG_', 'TTto2L2Nu', 'DYto2L', 'Zto2Nu', 'WtoLNu',
    'WW_', 'WZ_', 'ZZ_', 'Tbar', 'TW', 'ggH', 'VBF'
]

grouped_values = defaultdict(list)
others = {}

for s, vals in sample_values.items():
    matched = False
    for p in group_prefixes:
        if s.startswith(p):
            key = p.rstrip('_')
            grouped_values[key].extend(vals)
            matched = True
            break
    if not matched:
        others[s] = vals

# print grouped results first (average, min, max in minutes, and count)
print("# Grouped samples")
for key, vals in sorted(grouped_values.items(), key=lambda x: -sum(x[1])):
    count = len(vals)
    minutes_list = [v / 60.0 for v in vals]
    avg_minutes = sum(minutes_list) / count if count else 0.0
    sorted_mins = sorted(minutes_list) if count else []
    min_minutes = sorted_mins[0] if count else 0.0
    second_min_minutes = sorted_mins[1] if count > 1 else min_minutes
    max_minutes = sorted_mins[-1] if count else 0.0
    print(f"{key}\tavg={avg_minutes:.2f}\tmin={min_minutes:.2f}\t2nd_min={second_min_minutes:.2f}\tmax={max_minutes:.2f}\tn={count}")

print()
print("# Individual samples (not grouped)")
for s, vals in sorted(others.items(), key=lambda x: -sum(x[1])):
    count = len(vals)
    minutes_list = [v / 60.0 for v in vals]
    avg_minutes = sum(minutes_list) / count if count else 0.0
    sorted_mins = sorted(minutes_list) if count else []
    min_minutes = sorted_mins[0] if count else 0.0
    second_min_minutes = sorted_mins[1] if count > 1 else min_minutes
    max_minutes = sorted_mins[-1] if count else 0.0
    print(f"{s}\tavg={avg_minutes:.2f}\tmin={min_minutes:.2f}\t2nd_min={second_min_minutes:.2f}\tmax={max_minutes:.2f}\tn={count}")