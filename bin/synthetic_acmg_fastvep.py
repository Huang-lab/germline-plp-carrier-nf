#!/usr/bin/env python3
"""Test-only: synthesize an acmg_plp.tsv directly from a fixture VCF.

Stand-in for a real `fastvep annotate --acmg` + parse_fastvep_acmg.py run when
fastVEP and its supplementary databases are unavailable (the `test` profile,
gated behind params.allow_stub_acmg). Emits the SAME schema the real parser
does, so the DAG (carrier matrix, per-gene QC) is exercised unchanged:

  chr  pos  ref  alt  gene  acmg_label  acmg_criteria
  n_pathogenic_criteria  n_benign_criteria  is_acmg_PLP

Reads the inline SYN=... encoding written by bin/synthetic_vep.py (SYN subfield
8, if present, carries a synthetic ACMG label: 'Pathogenic'/'Likely pathogenic'/
'Benign'; absent => Uncertain significance), and falls back to the CSQ SYMBOL.
"""
from __future__ import annotations
import argparse
import gzip
import io
import sys

_HEADER = ("chr\tpos\tref\talt\tgene\tacmg_label\tacmg_criteria\t"
           "n_pathogenic_criteria\tn_benign_criteria\tis_acmg_PLP\n")

# synthetic label -> (criteria list) consistent with pathogenic/benign counts
_CRITERIA = {
    "pathogenic":        ["PVS1", "PM2_Supporting", "PP3"],
    "likely pathogenic": ["PM2_Supporting", "PP3"],
    "uncertain significance": [],
    "likely benign":     ["BP4"],
    "benign":            ["BA1"],
}
_LABEL_NORM = {
    "pathogenic": "Pathogenic",
    "likely pathogenic": "Likely_pathogenic",
    "uncertain significance": "Uncertain_significance",
    "likely benign": "Likely_benign",
    "benign": "Benign",
}
_PLP = {"pathogenic", "likely pathogenic"}
_PATHO = ("PVS", "PS", "PM", "PP")
_BENIGN = ("BA", "BS", "BP")


def _open(path: str):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with _open(args.inp) as fh, open(args.out, "w", encoding="utf-8") as out:
        out.write(_HEADER)
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt = parts[:5]
            info = parts[7]
            gene = ""
            label = "uncertain significance"
            for kv in info.split(";"):
                if kv.startswith("SYN="):
                    subs = kv[4:].split("|")
                    if subs and subs[0]:
                        gene = subs[0]
                    if len(subs) > 8 and subs[8]:
                        label = subs[8].strip().lower()
                    break
                if kv.startswith("CSQ="):
                    first = kv[4:].split(",", 1)[0].split("|")
                    if len(first) > 3:
                        gene = first[3]
            codes = _CRITERIA.get(label, [])
            npath = sum(1 for c in codes if c.startswith(_PATHO))
            nben = sum(1 for c in codes if c.startswith(_BENIGN))
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t"
                      f"{_LABEL_NORM.get(label, 'Uncertain_significance')}\t"
                      f"{';'.join(codes)}\t{npath}\t{nben}\t{int(label in _PLP)}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
