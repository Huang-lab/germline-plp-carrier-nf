#!/usr/bin/env python3
"""Parse a fastVEP `--acmg` VCF into the pipeline's acmg_plp.tsv schema.

fastVEP (Huang-lab/fastVEP) does its own consequence + supplementary annotation
+ ACMG-AMP classification and, in VCF output, appends two subfields to the CSQ
INFO field:
  * ACMG           — 5-tier shorthand: P / LP / VUS / LB / B
  * ACMG_CRITERIA  — met criteria codes joined by '&' (e.g. PVS1&PM2_Supporting)

This script collapses fastVEP's per-transcript CSQ blocks to one row per
variant (keeping the MOST SEVERE ACMG call across transcripts) and emits the
same columns the Nextflow pipeline's acmg_plp.tsv uses, so downstream
(carrier matrix, per-gene QC) consumes it unchanged:

  chr  pos  ref  alt  gene  acmg_label  acmg_criteria
  n_pathogenic_criteria  n_benign_criteria  is_acmg_PLP

Reuses plp_rules.csq for CSQ header parsing / value decoding.
"""
from __future__ import annotations
import argparse
import gzip
import io
import os as _os
import sys as _sys

_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.realpath(__file__))))
from plp_rules.csq import parse_csq_format, decode_csq_field  # noqa: E402

# 5-tier severity ordering (higher = more pathogenic) and shorthand→label map.
_SEVERITY = {"P": 4, "LP": 3, "VUS": 2, "LB": 1, "B": 0}
_LABEL = {
    "P": "Pathogenic",
    "LP": "Likely_pathogenic",
    "VUS": "Uncertain_significance",
    "LB": "Likely_benign",
    "B": "Benign",
}
_PLP = {"P", "LP"}
_PATHO_PREFIXES = ("PVS", "PS", "PM", "PP")
_BENIGN_PREFIXES = ("BA", "BS", "BP")


def _open(path: str):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _count_criteria(criteria_codes: list[str]) -> tuple[int, int]:
    npath = sum(1 for c in criteria_codes if c.startswith(_PATHO_PREFIXES))
    nben = sum(1 for c in criteria_codes if c.startswith(_BENIGN_PREFIXES))
    return npath, nben


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="fastVEP --acmg VCF output")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    schema = None
    i_sym = i_acmg = i_crit = -1

    def _opt(name):
        return schema.fields.index(name) if (schema and name in schema.fields) else -1

    with _open(args.vcf) as fh, open(args.out, "w", encoding="utf-8") as out:
        out.write("chr\tpos\tref\talt\tgene\tacmg_label\tacmg_criteria\t"
                  "n_pathogenic_criteria\tn_benign_criteria\tis_acmg_PLP\n")
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ,"):
                schema = parse_csq_format(line)
                i_sym = _opt("SYMBOL")
                i_acmg = _opt("ACMG")
                i_crit = _opt("ACMG_CRITERIA")
                if i_acmg < 0:
                    print("ERROR: no ACMG subfield in CSQ — was fastVEP run with --acmg?",
                          file=_sys.stderr)
                    return 2
                continue
            if line.startswith("#"):
                continue
            if schema is None:
                print("ERROR: no CSQ header before records", file=_sys.stderr)
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

            # Most-severe ACMG across transcript blocks.
            best_rank = -1
            best_sh = ""
            best_gene = ""
            best_crit = ""
            for block in csq.split(","):
                f = block.split("|")

                def g(idx):
                    return decode_csq_field(f[idx]) if 0 <= idx < len(f) else ""
                sh = g(i_acmg).strip()
                rank = _SEVERITY.get(sh, -1)
                if rank > best_rank:
                    best_rank = rank
                    best_sh = sh
                    best_gene = g(i_sym)
                    best_crit = g(i_crit)

            if best_rank < 0:
                continue  # no ACMG call on any transcript
            codes = [c for c in best_crit.replace("&", ";").split(";") if c]
            npath, nben = _count_criteria(codes)
            label = _LABEL.get(best_sh, best_sh)
            is_plp = int(best_sh in _PLP)
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{best_gene}\t{label}\t"
                      f"{';'.join(codes)}\t{npath}\t{nben}\t{is_plp}\n")
    return 0


if __name__ == "__main__":
    _sys.exit(main())
