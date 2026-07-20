#!/usr/bin/env python3
"""Classify AlphaMissense-calibrated P/LP.

Writes: chr, pos, ref, alt, gene, am_score, min_strength, is_AM_PLP.
"""
from __future__ import annotations
import argparse
import gzip
import io
import sys

from plp_rules.config import AlphaMissenseParams
from plp_rules.alphamissense import is_plp, load_calibration_tsv, normalize_strength
from plp_rules.csq import parse_csq_format, decode_csq_field


def _open(path: str):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--calibration", default="", help="Gene-specific AM threshold TSV")
    ap.add_argument("--min-strength", default="PP3_Moderate")
    ap.add_argument("--default-threshold", type=float, default=0.564)
    args = ap.parse_args()

    strength = normalize_strength(args.min_strength)
    calib = load_calibration_tsv(args.calibration) if args.calibration else {}
    params = AlphaMissenseParams(
        calibration_table=args.calibration,
        min_evidence_strength=strength,
        default_threshold=args.default_threshold,
    )
    schema = None
    idx_sym = idx_am = -1

    with _open(args.vcf) as fh, open(args.out, "w", encoding="utf-8") as out:
        out.write("chr\tpos\tref\talt\tgene\tam_score\tmin_strength\tis_AM_PLP\n")
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ,"):
                schema = parse_csq_format(line)
                idx_sym = schema.index_of("SYMBOL")
                idx_am = schema.index_of("am_pathogenicity")
                continue
            if line.startswith("#"):
                continue
            if schema is None:
                print("ERROR: no CSQ header before records", file=sys.stderr)
                return 2
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, _q, _f, info = parts[:8]
            csq = ""
            for kv in info.split(";"):
                if kv.startswith("CSQ="):
                    csq = kv[4:]
                    break
            if not csq:
                continue
            first = csq.split(",", 1)[0].split("|")
            def _g(i: int) -> str:
                return decode_csq_field(first[i]) if i < len(first) else ""
            gene = _g(idx_sym)
            am_raw = _g(idx_am)
            try:
                am = float(am_raw) if am_raw not in ("", ".", "NA") else None
            except ValueError:
                am = None
            plp = is_plp(gene, am, calib, params)
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{am_raw}\t{strength}\t{int(plp)}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
