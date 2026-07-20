#!/usr/bin/env python3
"""Validate a VEP-annotated VCF chunk: CSQ Format subfields + SYMBOL population.

Reads the input VCF (bgzipped ok — uses gzip.open transparently), extracts the
CSQ INFO header Description, and checks required subfields exist and SYMBOL is
populated in at least `min_symbol_fraction` of records.

Exit non-zero on failure; write a JSON report to stdout regardless.
"""
from __future__ import annotations
import argparse
import gzip
import io
import json
import sys
from typing import Iterator, Optional

from plp_rules.config import ValidateParams
from plp_rules.csq import validate_chunk


def _open(path: str):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _iter_vcf(path: str) -> Iterator[tuple[str, Optional[str]]]:
    """Yield ('#header', description_or_None) then ('record', csq_value or None) per line."""
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ,"):
                yield ("csq_header", line.rstrip("\n"))
                continue
            if line.startswith("#"):
                yield ("header", line.rstrip("\n"))
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            info = parts[7]
            csq: Optional[str] = None
            for kv in info.split(";"):
                if kv.startswith("CSQ="):
                    csq = kv[4:]
                    break
            yield ("record", csq)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--min-symbol-fraction", type=float, default=None)
    ap.add_argument("--required", nargs="*", default=None)
    args = ap.parse_args()

    csq_header_desc: Optional[str] = None
    csqs: list[str] = []
    for kind, val in _iter_vcf(args.vcf):
        if kind == "csq_header":
            csq_header_desc = val
        elif kind == "record":
            csqs.append(val or "")

    if csq_header_desc is None:
        print(json.dumps({"ok": False, "errors": ["No CSQ INFO header found"]}))
        return 2

    kw = {}
    if args.min_symbol_fraction is not None:
        kw["min_symbol_fraction"] = args.min_symbol_fraction
    if args.required:
        kw["required_csq_fields"] = tuple(args.required)
    params = ValidateParams(**kw)
    ok, errs = validate_chunk(csq_header_desc, csqs, params)
    print(json.dumps({"ok": ok, "errors": errs, "n_records": len(csqs)}))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
