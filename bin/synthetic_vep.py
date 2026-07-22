#!/usr/bin/env python3
"""Test-only: synthesize a VEP-like CSQ header + INFO annotations on a fixture VCF.

Reads the fixture VCF, expects each record's INFO to carry an inline `SYN=` field
with pipe-delimited synthetic annotations:
    SYN=<SYMBOL>|<Consequence>|<am_pathogenicity>|<LoF>|<CLNSIG>|<CLNREVSTAT>|<gnomAD_AF>|<gnomAD_AF_grpmax>
Rewrites the INFO to add a compatible CSQ tag so downstream modules can run unchanged.
"""
from __future__ import annotations
import argparse
import sys

CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
    'Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|MANE_SELECT|LoF|LoF_filter|'
    'am_pathogenicity|am_class|ClinVar_CLNSIG|ClinVar_CLNREVSTAT|ClinVar_CLNDN|ClinVar_CLNSIGCONF|'
    'gnomAD_AF|gnomAD_AF_grpmax">'
)


def build_csq(alt: str, syn: str) -> str:
    parts = syn.split("|")
    while len(parts) < 8:
        parts.append("")
    symbol, cons, am, lof, cnsig, crev, gaf, gmax = parts[:8]
    # VEP encodes commas in CSQ subfield values as '&' (its list separator).
    # Match that behavior in the fixture so decoders are exercised.
    def enc(s: str) -> str:
        return s.replace(",", "&")
    # Synthetic condition (CLNDN) present whenever there is a ClinVar sig;
    # CLNSIGCONF only when the call is Conflicting.
    clndn = "Synthetic_condition" if cnsig else ""
    clnsigconf = "Pathogenic(2)&Benign(1)" if cnsig and "Conflicting" in cnsig else ""
    fields = [alt, cons, "MODERATE", symbol, f"ENSG_{symbol}", "Transcript", f"ENST_{symbol}",
              "protein_coding", "1", lof, "", am, "",
              enc(cnsig), enc(crev), enc(clndn), clnsigconf, gaf, gmax]
    return "|".join(fields)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.inp, "r", encoding="utf-8") as fh, open(args.out, "w", encoding="utf-8") as out:
        emitted_csq_header = False
        for line in fh:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                if not emitted_csq_header:
                    out.write(CSQ_HEADER + "\n")
                    emitted_csq_header = True
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                out.write(line)
                continue
            alt = parts[4]
            info = parts[7]
            syn = ""
            new_kvs = list(info.split(";"))  # keep SYN so downstream fixtures (synthetic_intervar) can read the ACMG label
            for kv in info.split(";"):
                if kv.startswith("SYN="):
                    syn = kv[4:]
                    break
            csq = build_csq(alt, syn)
            new_kvs.append(f"CSQ={csq}")
            parts[7] = ";".join([kv for kv in new_kvs if kv])
            out.write("\t".join(parts) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
