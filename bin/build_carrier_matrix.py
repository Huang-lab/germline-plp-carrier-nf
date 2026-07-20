#!/usr/bin/env python3
"""Build the per-person carrier matrix (long TSV; optional wide pivot).

Inputs:
  --gt        bcftools query TSV: chr\tpos\tref\talt\tsample\tGT (one row per sample-variant carrying alt)
  --clinvar   classify_clinvar output TSV
  --acmg      acmg_postprocess output TSV
  --am        classify_alphamissense output TSV
  --keep      sample keep-list, one ID per line (optional)
  --out-long  long-format output TSV (chr,pos,ref,alt,gene,person_id,is_clinvar_PLP,is_acmg_PLP,is_AM_PLP)
  --out-wide  optional wide pivot (variant_key rows × person columns, cell = OR of the three flags)
"""
from __future__ import annotations
import argparse
import csv
from collections import defaultdict


def _key(chrom: str, pos: str, ref: str, alt: str) -> tuple[str, str, str, str]:
    return (chrom, pos, ref, alt)


def _load_classification(path: str, flag_col: str) -> tuple[dict, dict]:
    """Return {variant_key: 0/1 flag}, {variant_key: gene}."""
    flags: dict = {}
    genes: dict = {}
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            k = _key(row["chr"], row["pos"], row["ref"], row["alt"])
            flags[k] = int(row[flag_col])
            g = row.get("gene", "") or genes.get(k, "")
            if g:
                genes[k] = g
    return flags, genes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gt", required=True)
    ap.add_argument("--clinvar", default="")
    ap.add_argument("--acmg", default="")
    ap.add_argument("--am", default="")
    ap.add_argument("--keep", default="")
    ap.add_argument("--out-long", required=True)
    ap.add_argument("--out-wide", default="")
    args = ap.parse_args()

    cv_flags, cv_genes = (_load_classification(args.clinvar, "is_clinvar_PLP") if args.clinvar else ({}, {}))
    ac_flags, ac_genes = (_load_classification(args.acmg, "is_acmg_PLP") if args.acmg else ({}, {}))
    am_flags, am_genes = (_load_classification(args.am, "is_AM_PLP") if args.am else ({}, {}))

    keep: set[str] | None = None
    if args.keep:
        try:
            with open(args.keep, "r", encoding="utf-8") as fh:
                lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
            # Ignore obviously-non-ID single-line "placeholder"/"none" markers.
            if lines and lines != ["none"] and lines != ["all"]:
                keep = set(lines)
        except OSError:
            keep = None

    qualifying: set = set()
    for k, v in cv_flags.items():
        if v:
            qualifying.add(k)
    for k, v in ac_flags.items():
        if v:
            qualifying.add(k)
    for k, v in am_flags.items():
        if v:
            qualifying.add(k)

    def gene_for(k) -> str:
        return cv_genes.get(k) or ac_genes.get(k) or am_genes.get(k) or ""

    long_rows: list[list[str]] = []
    wide_cells: dict = defaultdict(dict)

    with open(args.gt, "r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 6:
                continue
            chrom, pos, ref, alt, sample, gt = row[:6]
            if keep is not None and sample not in keep:
                continue
            gt_norm = gt.replace("|", "/")
            if gt_norm in ("./.", ".", "0/0", ""):
                continue
            has_alt = any(a not in ("0", ".", "") for a in gt_norm.split("/"))
            if not has_alt:
                continue
            k = _key(chrom, pos, ref, alt)
            if k not in qualifying:
                continue
            cv = cv_flags.get(k, 0)
            ac = ac_flags.get(k, 0)
            am = am_flags.get(k, 0)
            gene = gene_for(k)
            long_rows.append([chrom, pos, ref, alt, gene, sample, str(cv), str(ac), str(am)])
            if args.out_wide:
                vk = f"{chrom}:{pos}:{ref}:{alt}"
                wide_cells[vk][sample] = 1 if (cv or ac or am) else wide_cells[vk].get(sample, 0)

    with open(args.out_long, "w", encoding="utf-8") as out:
        out.write("chr\tpos\tref\talt\tgene\tperson_id\tis_clinvar_PLP\tis_acmg_PLP\tis_AM_PLP\n")
        for r in long_rows:
            out.write("\t".join(r) + "\n")

    if args.out_wide:
        samples = sorted({s for cells in wide_cells.values() for s in cells})
        with open(args.out_wide, "w", encoding="utf-8") as out:
            out.write("variant_key\t" + "\t".join(samples) + "\n")
            for vk in sorted(wide_cells):
                cells = wide_cells[vk]
                out.write(vk + "\t" + "\t".join(str(cells.get(s, 0)) for s in samples) + "\n")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
