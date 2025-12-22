#!/usr/bin/env python3

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Optional


def _find_nominal_base(corrections_by_name: dict, explicit_nominal: Optional[str]) -> str:
    if explicit_nominal:
        if explicit_nominal not in corrections_by_name:
            raise KeyError(f"Nominal correction '{explicit_nominal}' not found in input JSON")
        base = explicit_nominal
        if base.endswith("_up") or base.endswith("_down"):
            raise ValueError("--nominal must be the base name (without _up/_down)")
        return base

    candidates: list[str] = []
    for name in corrections_by_name.keys():
        if name.endswith("_up") or name.endswith("_down"):
            continue
        if f"{name}_up" in corrections_by_name and f"{name}_down" in corrections_by_name:
            candidates.append(name)

    if not candidates:
        raise RuntimeError(
            "Could not auto-detect a nominal correction with matching _up/_down siblings. "
            "Pass --nominal <baseName>."
        )
    if len(candidates) > 1:
        raise RuntimeError(
            "Multiple nominal candidates found: "
            f"{', '.join(candidates)}. Pass --nominal <baseName>."
        )
    return candidates[0]


def convert_file(in_path: Path, out_path: Path, nominal_name: Optional[str], keep_originals: bool) -> None:
    doc = json.loads(in_path.read_text())
    if doc.get("schema_version") != 2:
        raise ValueError(f"Expected schema_version=2, got {doc.get('schema_version')!r}")

    corrections = doc.get("corrections")
    if not isinstance(corrections, list):
        raise ValueError("Input JSON does not contain a 'corrections' list")

    corrections_by_name = {c.get("name"): c for c in corrections}
    if None in corrections_by_name:
        raise ValueError("At least one correction is missing a 'name'")

    base = _find_nominal_base(corrections_by_name, nominal_name)
    corr_nom = corrections_by_name[base]
    corr_up = corrections_by_name[f"{base}_up"]
    corr_down = corrections_by_name[f"{base}_down"]

    # Build a single correction with an extra string input 'systematic'
    systematic_input = {
        "name": "systematic",
        "type": "string",
        "description": "Systematic variation (nominal, up, down)",
    }

    new_corr = {
        "name": corr_nom["name"],
        "description": (corr_nom.get("description") or "")
        + " (with systematic variations via 'systematic' input)",
        "version": corr_nom.get("version", 0),
        # Put 'systematic' first (common POG style), then keep original inputs.
        "inputs": [systematic_input, *copy.deepcopy(corr_nom.get("inputs", []))],
        "output": copy.deepcopy(corr_nom.get("output")),
        "data": {
            "nodetype": "category",
            "input": "systematic",
            "content": [
                {"key": "nominal", "value": copy.deepcopy(corr_nom.get("data"))},
                {"key": "up", "value": copy.deepcopy(corr_up.get("data"))},
                {"key": "down", "value": copy.deepcopy(corr_down.get("data"))},
            ],
        },
    }

    out_doc = dict(doc)
    if keep_originals:
        out_doc["corrections"] = [new_corr, corr_nom, corr_up, corr_down]
    else:
        out_doc["corrections"] = [new_corr]

    out_path.write_text(json.dumps(out_doc, indent=2, sort_keys=False) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Convert separate *_up/*_down corrections into a single CorrectionLib correction "
            "with a categorical 'systematic' input selecting nominal/up/down."
        )
    )
    ap.add_argument("input", type=Path, help="Input correctionlib JSON (schema_version=2)")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output JSON path (default: overwrite input when --in-place is set)",
    )
    ap.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite the input file (ignored if -o/--output is provided)",
    )
    ap.add_argument(
        "--nominal",
        default=None,
        help="Base correction name (without _up/_down). Needed if auto-detection is ambiguous.",
    )
    ap.add_argument(
        "--keep-originals",
        action="store_true",
        help="Keep the original three corrections in addition to the merged one.",
    )

    args = ap.parse_args()

    in_path: Path = args.input
    if not in_path.exists():
        raise FileNotFoundError(str(in_path))

    if args.output is not None:
        out_path = args.output
    else:
        if not args.in_place:
            ap.error("Provide -o/--output or use --in-place")
        out_path = in_path

    convert_file(in_path, out_path, args.nominal, args.keep_originals)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
