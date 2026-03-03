#!/usr/bin/env python3
"""Compute LHEScaleWeight inclusive sums for samples using dascache lists.

Reads dascache text files, opens a few NanoAOD root files per sample via xrootd,
computes the sum of the first 9 LHEScaleWeight entries, and writes a YAML
mapping usable as sampleCfg["LHEScaleSumw"].
"""

import argparse
import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Iterable, List, Optional, Tuple

import awkward as ak
import numpy as np
import uproot
import yaml

logger = logging.getLogger("compute_lhe_scale_sumw")

# Tame xrootd stalls for missing/offline files
os.environ.setdefault("XRD_REQUESTTIMEOUT", "20")
os.environ.setdefault("XRD_CONNECTTIMEOUT", "10")
os.environ.setdefault("XRD_CONNECTIONWINDOW", "15")
os.environ.setdefault("XRD_STREAMTIMEOUT", "20")


def to_file_path(path: str) -> Optional[str]:
    """Normalize a PFN to a readable path (prefer local /eos/cms)."""
    if not path:
        return None
    path = path.strip()
    if not path or path.startswith("#"):
        return None
    if path.startswith("/eos/cms"):
        path = path[len("/eos/cms") :]
    if path.startswith("/store"):
        return "/eos/cms" + path
    if path.startswith("root://cms-xrd-global.cern.ch/"):
        # Map AAA PFNs to the EOS mount if possible
        store_pos = path.find("/store")
        if store_pos != -1:
            return "/eos/cms" + path[store_pos:]
    if path.startswith("root://"):
        return path
    return None


def load_urls(list_path: str, max_files: int) -> List[str]:
    urls: List[str] = []
    with open(list_path) as f:
        for line in f:
            url = to_file_path(line)
            if not url:
                continue
            urls.append(url)
            if max_files > 0 and len(urls) >= max_files:
                break
    return urls


def sum_lhe_scale_weights(urls: Iterable[str], entry_stop: Optional[int]) -> Tuple[ak.Array, int]:
    """Sum first 9 LHEScaleWeight entries over provided files."""
    total = ak.Array(np.zeros(9, dtype="float64"))
    n_events = 0
    for url in urls:
        file_spec = None
        if url.startswith("/eos/") and not os.path.exists(url):
            # Use triple-slash after host to avoid "relative path disallowed"
            xrd_url = "root://cms-xrd-global.cern.ch///" + url[len("/eos/cms/") :].lstrip("/")
            # logger.warning("File not found locally, trying xrootd: %s", xrd_url)
            file_spec = {xrd_url: "Events"}
        else:
            file_spec = {url: "Events"}
        try:
            for arrays in uproot.iterate(
                file_spec,
                ["LHEScaleWeight"],
                step_size="100 MB",
                entry_stop=entry_stop,
            ):
                w = arrays["LHEScaleWeight"]
                n_events += len(w)
                w9 = ak.fill_none(ak.pad_none(w, 9)[:, :9], 0.0)
                total = total + ak.sum(w9, axis=0)
        except Exception:
            logger.exception("Failed to read %s", url)
    return total, n_events


def is_data_sample(sample: str) -> bool:
    low = sample.lower()
    return any(tok in low for tok in ["run202", "egamma", "muon", "jetmet", "single" ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dascache-dir",
        default=os.path.join(os.path.dirname(__file__), "..", "config", "dascache"),
        help="Directory with dascache txt files",
    )
    parser.add_argument(
        "--years",
        nargs="+",
        default=["2022", "2022EE", "2023", "2023BPix"],
        help="Years to include (matched inside filename)",
    )
    parser.add_argument(
        "--sample-filter",
        default=None,
        help="Optional substring to select samples (case-sensitive)",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=2,
        help="Max number of ROOT files per sample",
    )
    parser.add_argument(
        "--entry-stop",
        type=int,
        default=None,
        help="Optional entry_stop per file (default: all events)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel workers over samples (use >1 to speed up I/O)",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(os.path.dirname(__file__), "..", "data", "LHEScaleSumw.yaml"),
        help="Output YAML path",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    dascache_dir = os.path.abspath(args.dascache_dir)
    if not os.path.isdir(dascache_dir):
        raise SystemExit(f"dascache dir not found: {dascache_dir}")

    results = {}
    tasks = []
    files = sorted(os.listdir(dascache_dir))
    for fname in files:
        if not fname.endswith(".txt"):
            continue
        if args.sample_filter and args.sample_filter not in fname:
            continue
        if not any(year in fname for year in args.years):
            continue
        sample = os.path.splitext(fname)[0]
        if is_data_sample(sample):
            continue

        list_path = os.path.join(dascache_dir, fname)
        urls = load_urls(list_path, args.max_files)
        if not urls:
            logger.warning("No URLs found for %s", sample)
            continue

        tasks.append((sample, urls))

    if args.workers <= 1:
        for sample, urls in tasks:
            logger.info("Processing %s (%d files)", sample, len(urls))
            total, n_events = sum_lhe_scale_weights(urls, args.entry_stop)
            if n_events == 0:
                logger.warning("No events read for %s", sample)
                continue

            results[sample] = [float(x) for x in ak.to_numpy(total)]
            logger.info("  events: %d", n_events)
            logger.info("  sumw[:9]: %s", results[sample])
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_map = {
                executor.submit(sum_lhe_scale_weights, urls, args.entry_stop): (sample, urls)
                for sample, urls in tasks
            }

            for future in as_completed(future_map):
                sample, urls = future_map[future]
                try:
                    total, n_events = future.result()
                except Exception:
                    logger.exception("Sample %s failed", sample)
                    continue

                if n_events == 0:
                    logger.warning("No events read for %s", sample)
                    continue

                results[sample] = [float(x) for x in ak.to_numpy(total)]
                logger.info("Processing %s (%d files)", sample, len(urls))
                logger.info("  events: %d", n_events)
                logger.info("  sumw[:9]: %s", results[sample])

    if not results:
        logger.warning("No results produced; exiting without writing YAML")
        return

    out_path = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as out:
        yaml.safe_dump({"LHEScaleSumw": results}, out, sort_keys=True)
    logger.info("Wrote %d samples to %s", len(results), out_path)


if __name__ == "__main__":
    main()
