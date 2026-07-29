#!/usr/bin/env python3
"""Classify ClinVar P/LP from a VEP-annotated VCF chunk.

Writes a TSV: chr, pos, ref, alt, gene, clnsig, clnrevstat, condition,
clnsigconf, stars, is_clinvar_PLP.
"""
from __future__ import annotations
import argparse
import gzip
import io
import sys

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.realpath(__file__))))
from plp_rules.config import ClinVarParams
from plp_rules.clinvar import parse_stars, is_plp, parse_condition
from plp_rules.csq import parse_csq_format, decode_csq_field


def _open(path: str):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-stars", type=int, default=2)
    args = ap.parse_args()

    params = ClinVarParams(min_stars=args.min_stars)
    schema = None
    idx_sym = idx_sig = idx_rev = idx_dn = idx_conf = -1

    def _opt_idx(sch, name):
        # Tolerant: CLNDN/CLNSIGCONF may be absent in older ClinVar releases.
        return sch.fields.index(name) if name in sch.fields else -1

    with _open(args.vcf) as fh, open(args.out, "w", encoding="utf-8") as out:
        out.write("chr\tpos\tref\talt\tgene\tclnsig\tclnrevstat\tcondition\tclnsigconf\tstars\tis_clinvar_PLP\n")
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ,"):
                schema = parse_csq_format(line)
                idx_sym = schema.index_of("SYMBOL")
                idx_sig = schema.index_of("ClinVar_CLNSIG")
                idx_rev = schema.index_of("ClinVar_CLNREVSTAT")
                idx_dn = _opt_idx(schema, "ClinVar_CLNDN")
                idx_conf = _opt_idx(schema, "ClinVar_CLNSIGCONF")
                continue
            if line.startswith("#"):
                continue
            if schema is None:
                print("ERROR: no CSQ header before records", file=sys.stderr)
                return 2
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, _qual, _filt, info = parts[:8]
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
            clnsig = _g(idx_sig)
            revstat = _g(idx_rev)
            condition = parse_condition(_g(idx_dn)) if idx_dn >= 0 else ""
            clnsigconf = _g(idx_conf) if idx_conf >= 0 else ""
            stars = parse_stars(revstat)
            plp = is_plp(clnsig, revstat, params)
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{clnsig}\t{revstat}\t{condition}\t{clnsigconf}\t{stars}\t{int(plp)}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
