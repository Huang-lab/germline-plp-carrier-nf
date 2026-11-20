#!/usr/bin/env python3
"""Pull the carriers of a gene list out of the per-chromosome pick results
into one table, with each carrier's variant annotation joined in.

Reads (never writes) <pick-root>/<prefix><chr>/acmg_carriers.tsv (prefix
results-carriers-pick- by default; results-carriers-fastvep_pick- for the
fastVEP v0.5.0 fastvep_pick output)
and acmg_variants.tsv for every chromosome present, keeps rows whose gene is
in --genes, and writes:

  --out            one row per (carrier, variant, gene): the carrier columns,
                   then the variant table's annotation for that (variant, gene),
                   then n_plp_variants_in_gene = how many distinct P/LP variants
                   this person carries in this gene (2+ flags possible compound
                   heterozygotes, e.g. for a recessive gene such as MUTYH;
                   phase is unknown)
  --out-summary    per gene: P/LP variants, carriers, het/hom/hemi carriers,
                   carriers with 2+ P/LP variants in the gene

Gene names are matched exactly against the `gene` column (SYMBOL). Genes in
the list with no carrier, and genes not seen in any variant table at all
(possible typo or alias), are reported on stderr.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import carrier_markers  # noqa: E402

csv.field_size_limit(sys.maxsize)

CARRIER_KEY = ("chr", "pos", "ref", "alt", "gene")


def read_genes(path: str) -> list[str]:
    genes = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            for tok in line.replace(",", " ").split():
                if tok.startswith("#"):
                    break
                if tok not in genes:
                    genes.append(tok)
    return genes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pick-root", required=True,
                    help="dir holding results-carriers-pick-<chr>/ folders")
    ap.add_argument("--genes", required=True,
                    help="gene symbols, one per line (or space/comma separated); # comments ok")
    ap.add_argument("--prefix", default="results-carriers-pick-",
                    help="per-chromosome dir prefix (results-carriers-fastvep_pick- for "
                         "fastvep_pick output)")
    ap.add_argument("--allow-pick-anomalies", "--allow-failed", dest="allow_failed",
                    action="store_true",
                    help="read chromosome dirs holding a FAILED marker anyway (default: refuse)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args(argv)

    genes = read_genes(args.genes)
    gene_set = set(genes)
    if not genes:
        ap.error(f"no genes in {args.genes}")

    dirs = sorted(d for d in glob.glob(os.path.join(args.pick_root, args.prefix + "*"))
                  if os.path.isfile(os.path.join(d, "acmg_carriers.tsv")))
    if not dirs:
        print(f"ERROR: no {args.prefix}*/ under {args.pick_root}", file=sys.stderr)
        return 3
    for d in dirs:
        carrier_markers.check_dir(d, args.allow_failed, "extract_gene_carriers")
    chroms = [os.path.basename(d)[len(args.prefix):] for d in dirs]
    print(f"[extract] {len(genes)} gene(s); {len(dirs)} chromosome(s): {chroms}", file=sys.stderr)

    variant_cols: list[str] = []
    carrier_cols: list[str] = []
    seen_genes: set[str] = set()
    rows_out = []
    for d in dirs:
        # Variant annotation for this list's genes only (P/LP rows are the only
        # ones with carriers, but keep any row of a listed gene for the join).
        ann: dict[tuple, dict] = {}
        with open(os.path.join(d, "acmg_variants.tsv"), encoding="utf-8", newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if not variant_cols:
                variant_cols = [c for c in rd.fieldnames if c not in CARRIER_KEY]
            for r in rd:
                if r["gene"] in gene_set:
                    seen_genes.add(r["gene"])
                    ann[tuple(r[k] for k in CARRIER_KEY)] = r
        with open(os.path.join(d, "acmg_carriers.tsv"), encoding="utf-8", newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if not carrier_cols:
                carrier_cols = list(rd.fieldnames)
            for r in rd:
                if r["gene"] not in gene_set:
                    continue
                v = ann.get(tuple(r[k] for k in CARRIER_KEY))
                if v is None:
                    print(f"ERROR: carrier row with no variant row in {d}: "
                          f"{[r[k] for k in CARRIER_KEY]}", file=sys.stderr)
                    return 4
                rows_out.append((r, v))

    # distinct P/LP variants per (person, gene)
    per_person_gene = defaultdict(set)
    for r, _ in rows_out:
        per_person_gene[(r["person_id"], r["gene"])].add(tuple(r[k] for k in CARRIER_KEY[:4]))

    header = carrier_cols + variant_cols + ["n_plp_variants_in_gene"]
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r, v in rows_out:
            w.writerow([r[c] for c in carrier_cols] + [v[c] for c in variant_cols]
                       + [len(per_person_gene[(r["person_id"], r["gene"])])])

    summ = {g: {"variants": set(), "carriers": set(), "het": set(), "hom": set(), "hemi": set()}
            for g in genes}
    for r, _ in rows_out:
        s = summ[r["gene"]]
        s["variants"].add(tuple(r[k] for k in CARRIER_KEY[:4]))
        s["carriers"].add(r["person_id"])
        if r["zygosity"] in ("het", "hom", "hemi"):
            s[r["zygosity"]].add(r["person_id"])
    with open(args.out_summary, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "plp_variants_with_carriers", "carriers", "het_carriers",
                    "hom_carriers", "hemi_carriers", "carriers_2plus_plp_variants"])
        for g in genes:
            s = summ[g]
            multi = sum(1 for (p, gg), vs in per_person_gene.items() if gg == g and len(vs) > 1)
            w.writerow([g, len(s["variants"]), len(s["carriers"]), len(s["het"]),
                        len(s["hom"]), len(s["hemi"]), multi])

    none = [g for g in genes if not summ[g]["carriers"]]
    unseen = [g for g in genes if g not in seen_genes]
    persons = {r["person_id"] for r, _ in rows_out}
    print(f"[extract] {len(rows_out)} carrier rows, {len(persons)} unique carriers -> {args.out}",
          file=sys.stderr)
    if none:
        print(f"[extract] genes with no carriers: {none}", file=sys.stderr)
    if unseen:
        print(f"[extract] WARNING genes never seen in any variant table (typo/alias, or their "
              f"chromosome is not in {args.pick_root} yet): {unseen}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
