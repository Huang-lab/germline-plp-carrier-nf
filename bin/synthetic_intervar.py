#!/usr/bin/env python3
"""Test-only: synthesize an InterVar-shaped table from a VEP-annotated fixture VCF.

Uses the inline SYN=... encoding written by synthetic_vep.py — the last two
optional SYN subfields (positions 8, 9) may carry a synthetic ACMG label if
present, e.g. 'Pathogenic' or 'Likely pathogenic'. Absent => VUS.
"""
from __future__ import annotations
import argparse
import sys


HEADER = ["#Chr", "Start", "End", "Ref", "Alt", "Ref.Gene", "InterVar: InterVar and Evidence"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.inp, "r", encoding="utf-8") as fh, open(args.out, "w", encoding="utf-8") as out:
        out.write("\t".join(HEADER) + "\n")
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt = parts[:5]
            info = parts[7]
            symbol = ""
            acmg_label = "Uncertain significance"
            for kv in info.split(";"):
                if kv.startswith("SYN="):
                    subs = kv[4:].split("|")
                    if len(subs) > 0:
                        symbol = subs[0]
                    if len(subs) > 8 and subs[8]:
                        acmg_label = subs[8]
                    break
                if kv.startswith("CSQ="):
                    first = kv[4:].split(",", 1)[0].split("|")
                    if len(first) > 3:
                        symbol = first[3]
            iv = f"InterVar: {acmg_label} PVS1=0 PS=[0,0,0,0,0] PM=[0,0,0,0,0,0,0] PP=[0,0,0,0,0] BA1=0 BS=[0,0,0,0] BP=[0,0,0,0,0,0,0,0]"
            end = str(int(pos) + max(len(ref), len(alt)) - 1)
            out.write("\t".join([chrom, pos, end, ref, alt, symbol, iv]) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
