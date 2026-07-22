#!/usr/bin/env python3
"""Extract per-sample GTs from a (optionally bgzipped) VCF without bcftools.

Used by CARRIER_MATRIX as the production GT-extraction path so the process
needs only python3 (one container), not bcftools+python together. Handles
multi-sample real pVCF chunks.

Emits chr\tpos\tref\talt\tsample\tGT rows for records whose (chr,pos) fall in
the qualifying BED. Genotype QC (DP/GQ/AB) is already applied upstream in
NORM_QC, so only GT is read here."""
from __future__ import annotations
import argparse
import gzip
import io
import sys


def _open(path: str):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _load_positions(path: str) -> set[tuple[str, int]]:
    """Load a CHROM<TAB>POS (1-based) positions file — same format bcftools -T uses."""
    keep = set()
    if not path:
        return keep
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            keep.add((cols[0], int(cols[1])))
    return keep


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--positions", required=True,
                    help="CHROM<TAB>POS (1-based) file, as used by bcftools -T")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    keep = _load_positions(args.positions)
    samples: list[str] = []

    with _open(args.vcf) as fh, open(args.out, "w", encoding="utf-8") as out:
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
