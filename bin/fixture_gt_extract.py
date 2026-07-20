#!/usr/bin/env python3
"""Test-only: extract per-sample GTs directly from a fixture VCF (no bcftools).

Emits chr\tpos\tref\talt\tsample\tGT rows for records whose (chr,pos) fall in
the qualifying BED. Ignores DP/GQ/AD (already applied earlier in the fallback
NORM_QC path)."""
from __future__ import annotations
import argparse
import sys


def _load_bed(path: str) -> set[tuple[str, int]]:
    keep = set()
    if not path:
        return keep
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c, s, e = line.rstrip("\n").split("\t")[:3]
            for p in range(int(s) + 1, int(e) + 1):
                keep.add((c, p))
    return keep


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    keep = _load_bed(args.bed)
    samples: list[str] = []

    with open(args.vcf, "r", encoding="utf-8") as fh, open(args.out, "w", encoding="utf-8") as out:
        for line in fh:
            if line.startswith("##"):
                continue
            parts = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = parts[9:]
                continue
            if len(parts) < 10:
                continue
            chrom, pos, _id, ref, alt = parts[:5]
            fmt = parts[8].split(":")
            try:
                gt_idx = fmt.index("GT")
            except ValueError:
                continue
            if (chrom, int(pos)) not in keep:
                continue
            for smpl, cell in zip(samples, parts[9:]):
                gt = cell.split(":")[gt_idx] if cell else "./."
                out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{smpl}\t{gt}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
