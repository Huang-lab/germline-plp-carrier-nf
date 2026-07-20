#!/usr/bin/env python3
"""Genotype QC helper. Reads a bcftools-query stream

    chr\tpos\tref\talt\tsample\tGT\tDP\tGQ\tAD

writes the same rows with GT set to ./. when QC thresholds fail.

The Nextflow norm_qc process typically prefers `bcftools +setGT` directly for
performance; this script exists so the QC rule (implemented once in
plp_rules.qc.genotype_passes) can also drive a fallback path or be tested.
"""
from __future__ import annotations
import argparse
import sys
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.realpath(__file__))))
from plp_rules.config import QCThresholds
from plp_rules.qc import genotype_passes


def _split_ad(ad: str) -> tuple[int, int]:
    if not ad or ad in (".", ""):
        return 0, 0
    parts = [int(x) if x not in ("", ".") else 0 for x in ad.split(",")]
    ref = parts[0] if parts else 0
    alt = parts[1] if len(parts) > 1 else 0
    return ref, alt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--min-gq", type=int, default=20)
    ap.add_argument("--het-ab-min", type=float, default=0.20)
    ap.add_argument("--het-ab-max", type=float, default=0.80)
    ap.add_argument("--hom-ab-min", type=float, default=0.90)
    args = ap.parse_args()

    th = QCThresholds(min_dp=args.min_dp, min_gq=args.min_gq,
                      het_ab_min=args.het_ab_min, het_ab_max=args.het_ab_max,
                      hom_ab_min=args.hom_ab_min)

    for line in sys.stdin:
        row = line.rstrip("\n").split("\t")
        if len(row) < 6:
            continue
        chrom, pos, ref, alt, sample, gt = row[:6]
        dp = int(row[6]) if len(row) > 6 and row[6] not in ("", ".") else None
        gq = int(row[7]) if len(row) > 7 and row[7] not in ("", ".") else None
        ad = row[8] if len(row) > 8 else ""
        ad_ref, ad_alt = _split_ad(ad)
        total = ad_ref + ad_alt
        ab = (ad_alt / total) if total > 0 else None
        gt_norm = gt.replace("|", "/")
        is_het = gt_norm in ("0/1", "1/0")
        is_hom_alt = gt_norm == "1/1"
        keep = True
        if is_het or is_hom_alt:
            keep = genotype_passes(dp, gq, ab, is_het, is_hom_alt, th)
        new_gt = gt if keep else "./."
        sys.stdout.write("\t".join([chrom, pos, ref, alt, sample, new_gt]) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
