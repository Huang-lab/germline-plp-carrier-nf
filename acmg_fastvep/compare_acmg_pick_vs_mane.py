#!/usr/bin/env python3
"""Compare MANE-only ACMG carrier results with per-gene pick results.

Reads (never writes) both result trees:
  --mane-root  dir holding results-carriers-<chr>/       (build_acmg_carriers.py)
  --pick-root  dir holding results-carriers-pick-<chr>/  (build_acmg_carriers_pick.py)
and writes TSVs into --out (normally <pick-root>/comparison/):

  overall.tsv               P/LP (variant, gene) rows, unique P/LP variants, carrier
                            rows, unique carriers (+ hemizygous), per method and delta
  per_gene.tsv              the same per gene (genes with P/LP in either method)
  per_chromosome.tsv        the same per chromosome
  pick_tier_counts.tsv      pick rows by pick_tier, all rows and P/LP rows
  mane_discordant_plp.tsv   pick rows with mane_discordant=1 where either call is P/LP
  plp_only_in_one.tsv       P/LP (variant, gene) rows found by one method only, with
                            the pick-side transcript/tier/is_mane_select for context

Only chromosomes present in both trees are compared; the rest are listed on
stderr. Streams the tables and keeps only P/LP rows in memory.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

PLP_LABELS = {"Pathogenic", "Likely_pathogenic"}


def _chr_dirs(root: str, prefix: str) -> dict[str, str]:
    out = {}
    for d in sorted(glob.glob(os.path.join(root, prefix + "*"))):
        name = os.path.basename(d.rstrip("/"))[len(prefix):]
        # "results-carriers-" also prefixes "results-carriers-pick-"; keep them apart.
        if prefix == "results-carriers-" and name.startswith("pick-"):
            continue
        if os.path.isfile(os.path.join(d, "acmg_variants.tsv")):
            out[name] = d
    return out


def _rows(path):
    with open(path, encoding="utf-8", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


class Side:
    def __init__(self):
        self.plp = {}                         # (chr,pos,ref,alt,gene) -> row
        self.carriers = defaultdict(set)      # gene -> {person}
        self.carrier_rows = defaultdict(int)  # gene -> n
        self.hemi = defaultdict(set)          # gene -> {person}
        self.by_chr = defaultdict(lambda: {"plp_rows": 0, "plp_variants": set(),
                                           "carriers": set(), "carrier_rows": 0})


def load(d: str, chrom: str, side: Side, tiers=None, discordant=None):
    for r in _rows(os.path.join(d, "acmg_variants.tsv")):
        if tiers is not None:
            t = r.get("pick_tier", "")
            plp = r["is_acmg_PLP"] == "1"
            tiers[t or "(sentinel)"][0] += 1
            tiers[t or "(sentinel)"][1] += int(plp)
            if r.get("mane_discordant") == "1" and (plp or r.get("acmg_label_mane") in PLP_LABELS):
                discordant.append(r)
        if r["is_acmg_PLP"] != "1":
            continue
        key = (r["chr"], r["pos"], r["ref"], r["alt"], r["gene"])
        side.plp[key] = r
        c = side.by_chr[chrom]
        c["plp_rows"] += 1
        c["plp_variants"].add(key[:4])
    for r in _rows(os.path.join(d, "acmg_carriers.tsv")):
        g, p = r["gene"], r["person_id"]
        side.carriers[g].add(p)
        side.carrier_rows[g] += 1
        if r.get("zygosity") == "hemi":
            side.hemi[g].add(p)
        c = side.by_chr[chrom]
        c["carriers"].add(p)
        c["carrier_rows"] += 1


def _write(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mane-root", required=True)
    ap.add_argument("--pick-root", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    mane_real = os.path.realpath(args.mane_root)
    out_real = os.path.realpath(args.out)
    if out_real == mane_real or out_real.startswith(mane_real + os.sep):
        ap.error("--out must not be inside --mane-root")

    mane_dirs = _chr_dirs(args.mane_root, "results-carriers-")
    pick_dirs = _chr_dirs(args.pick_root, "results-carriers-pick-")
    common = sorted(set(mane_dirs) & set(pick_dirs))
    for name, only in (("MANE", set(mane_dirs) - set(pick_dirs)),
                       ("pick", set(pick_dirs) - set(mane_dirs))):
        if only:
            print(f"[compare] only in {name}, NOT compared: {sorted(only)}", file=sys.stderr)
    if not common:
        print("[compare] ERROR: no chromosome present in both trees", file=sys.stderr)
        return 3
    print(f"[compare] comparing {len(common)} chromosome(s): {common}", file=sys.stderr)

    mane, pick = Side(), Side()
    tiers = defaultdict(lambda: [0, 0])
    discordant: list[dict] = []
    for chrom in common:
        load(mane_dirs[chrom], chrom, mane)
        load(pick_dirs[chrom], chrom, pick, tiers, discordant)

    os.makedirs(args.out, exist_ok=True)

    def totals(s: Side):
        persons = set().union(*s.carriers.values()) if s.carriers else set()
        hemi = set().union(*s.hemi.values()) if s.hemi else set()
        return [len(s.plp), len({k[:4] for k in s.plp}), sum(s.carrier_rows.values()),
                len(persons), len(hemi)]

    metrics = ["plp_variant_gene_rows", "plp_unique_variants", "carrier_rows",
               "unique_carriers", "unique_hemizygous_carriers"]
    tm, tp = totals(mane), totals(pick)
    _write(os.path.join(args.out, "overall.tsv"), ["metric", "mane_only", "pick", "delta"],
           [[m, a, b, b - a] for m, a, b in zip(metrics, tm, tp)])

    def plp_vars_by_gene(s: Side):
        d = defaultdict(set)
        for k in s.plp:
            d[k[4]].add(k[:4])
        return d

    mv, pv = plp_vars_by_gene(mane), plp_vars_by_gene(pick)
    genes = sorted(set(mv) | set(pv) | set(mane.carriers) | set(pick.carriers))
    rows = []
    for g in genes:
        a, b = len(mv.get(g, ())), len(pv.get(g, ()))
        ca, cb = len(mane.carriers.get(g, ())), len(pick.carriers.get(g, ()))
        rows.append([g, a, b, b - a, ca, cb, cb - ca, len(pick.hemi.get(g, ()))])
    rows.sort(key=lambda r: (-abs(r[3]), -abs(r[6]), r[0]))
    _write(os.path.join(args.out, "per_gene.tsv"),
           ["gene", "plp_variants_mane", "plp_variants_pick", "plp_variants_delta",
            "carriers_mane", "carriers_pick", "carriers_delta", "hemi_carriers_pick"], rows)

    rows = []
    for chrom in common:
        a, b = mane.by_chr[chrom], pick.by_chr[chrom]
        rows.append([chrom, a["plp_rows"], b["plp_rows"], len(a["plp_variants"]),
                     len(b["plp_variants"]), len(a["carriers"]), len(b["carriers"]),
                     a["carrier_rows"], b["carrier_rows"]])
    _write(os.path.join(args.out, "per_chromosome.tsv"),
           ["chromosome", "plp_rows_mane", "plp_rows_pick", "plp_variants_mane",
            "plp_variants_pick", "carriers_mane", "carriers_pick", "carrier_rows_mane",
            "carrier_rows_pick"], rows)

    total_rows = sum(v[0] for v in tiers.values()) or 1
    _write(os.path.join(args.out, "pick_tier_counts.tsv"),
           ["pick_tier", "rows", "fraction_of_rows", "plp_rows"],
           [[t, n, f"{n / total_rows:.6f}", p]
            for t, (n, p) in sorted(tiers.items(), key=lambda kv: -kv[1][0])])

    dcols = ["chr", "pos", "ref", "alt", "gene", "transcript_id", "pick_tier", "acmg_label",
             "acmg_label_mane", "acmg_criteria", "consequence", "hgvsc", "hgvsp"]
    _write(os.path.join(args.out, "mane_discordant_plp.tsv"), dcols,
           [[r.get(c, "") for c in dcols] for r in discordant])

    rows = []
    for key in sorted(set(mane.plp) ^ set(pick.plp)):
        in_mane = key in mane.plp
        r = mane.plp.get(key) or pick.plp[key]
        pr = pick.plp.get(key, {})
        rows.append([*key, "mane_only" if in_mane else "pick_only", r["acmg_label"],
                     r.get("consequence", ""), pr.get("transcript_id", ""),
                     pr.get("pick_tier", ""), pr.get("is_mane_select", "")])
    _write(os.path.join(args.out, "plp_only_in_one.tsv"),
           ["chr", "pos", "ref", "alt", "gene", "found_in", "acmg_label", "consequence",
            "pick_transcript_id", "pick_tier", "pick_is_mane_select"], rows)

    print("[compare] overall:", file=sys.stderr)
    for m, a, b in zip(metrics, tm, tp):
        print(f"  {m:28s} mane={a:>10} pick={b:>10} delta={b - a:+}", file=sys.stderr)
    print(f"  mane_discordant P/LP rows: {len(discordant)}; "
          f"P/LP rows in one method only: {len(rows)}", file=sys.stderr)
    print(f"[compare] wrote {args.out}/", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
