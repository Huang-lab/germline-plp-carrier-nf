#!/usr/bin/env python3
"""Post-process InterVar output → ACMG P/LP table.

Input: InterVar `*.hg38_multianno.txt.intervar` file (or equivalent) with a
column literally named `InterVar: InterVar and Evidence` or similar containing
the parseable label. We locate any column whose value starts with `InterVar:`
per row, so header wording variations are tolerated.

Emits: chr, pos, ref, alt, gene, acmg_label, acmg_criteria,
       n_pathogenic_criteria, n_benign_criteria, is_acmg_PLP.
"""
from __future__ import annotations
import argparse
import csv
import sys

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.realpath(__file__))))
from plp_rules.acmg import parse_intervar_evidence, reclassify


def _find_col(header: list[str], *candidates: str) -> int:
    lut = {h.strip().lower(): i for i, h in enumerate(header)}
    for c in candidates:
        i = lut.get(c.lower())
        if i is not None:
            return i
    return -1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--intervar", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--demote-on-ba1", action="store_true",
                    help="Force non-P/LP when the stand-alone benign criterion BA1 fired.")
    args = ap.parse_args()

    with open(args.intervar, "r", encoding="utf-8") as fh, open(args.out, "w", encoding="utf-8") as out:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        c_chr = _find_col(header, "#Chr", "Chr", "chr")
        c_start = _find_col(header, "Start")
        c_ref = _find_col(header, "Ref")
        c_alt = _find_col(header, "Alt")
        c_gene = _find_col(header, "Ref.Gene", "Gene.refGene", "gene")
        # InterVar's label column name has varied ("InterVar: InterVar and Evidence"); locate by value on first data row instead.
        c_iv = -1
        out.write("chr\tpos\tref\talt\tgene\tacmg_label\tacmg_criteria\t"
                  "n_pathogenic_criteria\tn_benign_criteria\tis_acmg_PLP\n")
        for row in reader:
            if c_iv < 0:
                for i, v in enumerate(row):
                    if isinstance(v, str) and v.startswith("InterVar:"):
                        c_iv = i
                        break
                if c_iv < 0:
                    continue
            ev = parse_intervar_evidence(row[c_iv] if c_iv < len(row) else "")
            if args.demote_on_ba1:
                ev = reclassify(ev, demote_on_benign_standalone=True)
            chrom = row[c_chr] if c_chr >= 0 else ""
            pos = row[c_start] if c_start >= 0 else ""
            ref = row[c_ref] if c_ref >= 0 else ""
            alt = row[c_alt] if c_alt >= 0 else ""
            gene = row[c_gene] if c_gene >= 0 else ""
            out.write(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{ev.label}\t{ev.criteria_str}\t"
                f"{len(ev.pathogenic_criteria)}\t{len(ev.benign_criteria)}\t{int(ev.is_plp)}\n"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
