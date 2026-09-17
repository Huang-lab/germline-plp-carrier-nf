#!/usr/bin/env python3
"""Retrieve fastVEP's own ACMG judgement + functional annotation per variant,
MANE Select transcript only, and match it to that same VCF's own genotype
columns to find carriers. Standalone reprocessing of already-annotated
`*.fastvep.vcf.gz` files (Huang-lab/fastVEP) — no fastVEP re-run.

Deliberately retrieval-only: no cross-transcript scoring or fallback chain.
For each variant, this looks up every CSQ block whose MANE_SELECT field is
set and reports exactly that block's already-computed ACMG call and
functional fields — **one output row per (variant, gene)**, not one row per
variant. At gene-dense loci a single position can overlap two or more genes,
each with its own independently MANE-tagged transcript (measured on a real
chr17 batch: 27.5% of variants there overlap 2+ genes this way); picking only
one of them arbitrarily would silently misattribute results exactly like the
original max-over-transcripts bug, just for a different reason. A variant
with no MANE-tagged transcript in its CSQ at all is NOT silently resolved
some other way — it is flagged (acmg_label=NO_MANE_TRANSCRIPT) so a human
decides what to do with it, matching neither fastVEP's own `--pick` fallback
chain nor the old max-over-transcripts behavior this replaces.

This exactly reproduces the highest-priority tier of fastVEP's own `--pick`
default order (`mane_select` first, see Huang-lab/fastVEP
crates/fastvep-annotate/src/pick.rs), so it reports what fastVEP's own
reduction mechanism would choose, without needing to re-run annotation with
--pick.

Carrier identification reads genotypes directly from the same
`.fastvep.vcf.gz`'s own FORMAT/sample columns (fastVEP passes them through
untouched) — no separate genotype VCF is needed. Sample columns are only
parsed for variants already classified P/LP, since that's a small fraction
of all variants and the cohort has tens of thousands of sample columns per
line.

Output keyed on chr/pos/ref/alt, matching bin/build_carrier_matrix.py's
convention, so a ClinVar table can be joined in later without reshaping.
"""
from __future__ import annotations
import argparse
import gzip
import io
import os as _os
import sys as _sys

_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.realpath(__file__))))
from plp_rules.csq import parse_csq_format, decode_csq_field  # noqa: E402
# Reuse the pipeline's existing ACMG label/severity tables and criteria
# counter rather than redefining them — one source of truth for what counts
# as a pathogenic/benign criterion code.
from parse_fastvep_acmg import _LABEL, _PLP, _count_criteria  # noqa: E402

NO_MANE_SENTINEL = "NO_MANE_TRANSCRIPT"

CSQ_FUNCTIONAL_FIELDS = ("Consequence", "IMPACT", "HGVSc", "HGVSp", "EXON", "INTRON")
SPLICEAI_FIELDS = ("DS_AG", "DS_AL", "DS_DG", "DS_DL")

# The standard fixed set of "predicted loss-of-function" consequence terms
# (gnomAD/LOFTEE convention: stop_gained, frameshift_variant, splice acceptor/
# donor, start_lost). This is a category-membership lookup against a term list
# used identically across the field, not a judgement call this pipeline makes
# — it does NOT account for NMD escape or last-exon position, so it is the
# basic pLoF flag, not the more rigorous LOFTEE-style one. See docs/ACMG_CRITERIA.md.
LOF_CONSEQUENCES = frozenset({
    "stop_gained", "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant",
    "start_lost",
})

VARIANT_HEADER = (
    "chr\tpos\tref\talt\tgene\tacmg_label\tacmg_criteria\t"
    "n_pathogenic_criteria\tn_benign_criteria\tis_acmg_PLP\t"
    "consequence\timpact\thgvsc\thgvsp\texon\tintron\tis_lof_consequence\t"
    "spliceai_ds_ag\tspliceai_ds_al\tspliceai_ds_dg\tspliceai_ds_dl\n"
)
CARRIER_HEADER = "chr\tpos\tref\talt\tgene\tperson_id\tGT\tzygosity\n"


def _open(path: str):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _find_mane_blocks_by_gene(blocks, i_mane, i_sym):
    """CSQ blocks whose MANE_SELECT field is non-empty, grouped by gene symbol.

    No scoring, no fallback: within one gene, MANE Select names exactly one
    transcript, so more than one MANE-tagged block for the same gene is a
    genuine anomaly (reported by the caller), not a tie to break. Across
    genes, an overlapping locus can and does legitimately have several genes
    each with their own MANE-tagged block — every one of them is returned,
    keyed by gene, so the caller reports each gene independently rather than
    picking one.
    """
    if i_mane < 0:
        return {}
    by_gene: dict[str, list[list[str]]] = {}
    for block in blocks:
        f = block.split("|")
        val = decode_csq_field(f[i_mane]).strip() if 0 <= i_mane < len(f) else ""
        if not val:
            continue
        gene = decode_csq_field(f[i_sym]).strip() if 0 <= i_sym < len(f) else ""
        by_gene.setdefault(gene, []).append(f)
    return by_gene


def _field(f, idx):
    return decode_csq_field(f[idx]).strip() if 0 <= idx < len(f) else ""


def _is_lof(consequence: str) -> int:
    """1 if any '&'-joined Consequence term is in the standard pLoF set."""
    return int(any(term in LOF_CONSEQUENCES for term in consequence.split("&")))


def _summarize_unmapped_blocks(blocks, i_sym, i_func_consequence):
    """For a variant with no MANE-tagged transcript: the union of gene symbols
    and consequence terms seen across every block, so a flagged row still says
    what's actually there (e.g. genuinely intergenic vs. a real gene lacking
    MANE) instead of going blank. Purely descriptive — no block is picked."""
    genes: list[str] = []
    consequences: list[str] = []
    for block in blocks:
        f = block.split("|")
        gene = _field(f, i_sym)
        if gene and gene not in genes:
            genes.append(gene)
        cons = _field(f, i_func_consequence)
        if cons and cons not in consequences:
            consequences.append(cons)
    return ";".join(genes), ";".join(consequences)


def _parse_spliceai_info(info: str, spliceai_id: str | None):
    """Return {allele: (DS_AG, DS_AL, DS_DG, DS_DL)} from the SpliceAI INFO
    field (format ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL),
    or {} if absent/not requested."""
    if not spliceai_id:
        return {}
    out = {}
    for kv in info.split(";"):
        if not kv.startswith(spliceai_id + "="):
            continue
        val = kv[len(spliceai_id) + 1:]
        for block in val.split(","):
            f = block.split("|")
            if len(f) < 6:
                continue
            allele = decode_csq_field(f[0])
            out[allele] = (f[2], f[3], f[4], f[5])
    return out


def _zygosity(gt: str):
    """(is_carrier, zygosity_label) for a VCF GT string. None for 0/0 or missing."""
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or "." in alleles:
        return None
    a, b = alleles
    if a == "0" and b == "0":
        return None
    return "hom" if a == b else "het"


def process_vcf(path: str, out_variants, out_carriers, stats: dict) -> None:
    schema = None
    i_sym = i_acmg = i_crit = i_mane = -1
    i_func = {}
    spliceai_id = None
    samples: list[str] = []
    gt_idx_cache = {}

    with _open(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=SpliceAI,"):
                spliceai_id = "SpliceAI"
                continue
            if line.startswith("##INFO=<ID=CSQ,"):
                schema = parse_csq_format(line)

                def _opt(name):
                    return schema.fields.index(name) if name in schema.fields else -1

                i_sym = _opt("SYMBOL")
                i_acmg = _opt("ACMG")
                i_crit = _opt("ACMG_CRITERIA")
                i_mane = _opt("MANE_SELECT")
                for name in CSQ_FUNCTIONAL_FIELDS:
                    i_func[name] = _opt(name)
                if i_acmg < 0:
                    print(f"ERROR: no ACMG subfield in CSQ ({path}) — was fastVEP run with --acmg?",
                          file=_sys.stderr)
                    _sys.exit(2)
                if i_mane < 0:
                    print(f"ERROR: no MANE_SELECT subfield in CSQ ({path}) — cannot do "
                          "MANE-only retrieval on this fastVEP build/CSQ schema", file=_sys.stderr)
                    _sys.exit(2)
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:] if len(cols) > 9 else []
                continue
            if line.startswith("#"):
                continue
            if schema is None:
                print(f"ERROR: no CSQ header before records ({path})", file=_sys.stderr)
                _sys.exit(2)

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

            stats["total_variants"] += 1
            blocks = csq.split(",")
            by_gene = _find_mane_blocks_by_gene(blocks, i_mane, i_sym)

            if not by_gene:
                stats["no_mane"] += 1
                genes_seen, consequences_seen = _summarize_unmapped_blocks(
                    blocks, i_sym, i_func["Consequence"])
                out_variants.write(
                    f"{chrom}\t{pos}\t{ref}\t{alt}\t{genes_seen}\t{NO_MANE_SENTINEL}\t\t0\t0\t0\t"
                    f"{consequences_seen}\t\t\t\t\t\t\t\t\t\t\n"
                )
                continue

            if len(by_gene) > 1:
                stats["multi_gene_variants"] += 1

            any_plp_gene = False
            plp_genes_this_variant = []
            for gene, gene_blocks in by_gene.items():
                stats["gene_rows"] += 1
                if len(gene_blocks) > 1:
                    # Same gene, more than one MANE_SELECT-tagged transcript: MANE
                    # names exactly one transcript per gene, so this is a genuine
                    # data anomaly, not a tie to break by picking one.
                    stats["same_gene_multi_mane_anomaly"] += 1
                    print(f"WARN: {chrom}:{pos} {ref}>{alt} gene={gene!r} has "
                          f"{len(gene_blocks)} MANE_SELECT-tagged transcripts for the "
                          "SAME gene; using the first by transcript order for this gene",
                          file=_sys.stderr)
                f = gene_blocks[0]

                sh = _field(f, i_acmg)
                crit = _field(f, i_crit)
                codes = [c for c in crit.replace("&", ";").split(";") if c]
                npath, nben = _count_criteria(codes)
                label = _LABEL.get(sh, sh)
                is_plp = int(sh in _PLP)
                func = {name: _field(f, i_func[name]) for name in CSQ_FUNCTIONAL_FIELDS}
                is_lof = _is_lof(func["Consequence"])

                spliceai_by_allele = _parse_spliceai_info(info, spliceai_id)
                ds = spliceai_by_allele.get(alt, ("", "", "", ""))

                out_variants.write(
                    f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{label}\t{';'.join(codes)}\t"
                    f"{npath}\t{nben}\t{is_plp}\t"
                    f"{func['Consequence']}\t{func['IMPACT']}\t{func['HGVSc']}\t{func['HGVSp']}\t"
                    f"{func['EXON']}\t{func['INTRON']}\t{is_lof}\t"
                    f"{ds[0]}\t{ds[1]}\t{ds[2]}\t{ds[3]}\n"
                )

                if is_plp == 1:
                    any_plp_gene = True
                    plp_genes_this_variant.append(gene)

            if not any_plp_gene or not samples:
                continue

            # Only variants with at least one P/LP gene call pay for
            # sample-column parsing, and a carrier row is written per
            # qualifying gene so a carrier of a variant flagged P/LP for two
            # different overlapping genes is recorded against both.
            format_field = parts[8] if len(parts) > 8 else ""
            if format_field not in gt_idx_cache:
                gt_fields = format_field.split(":")
                gt_idx_cache[format_field] = gt_fields.index("GT") if "GT" in gt_fields else -1
            gt_idx = gt_idx_cache[format_field]
            if gt_idx < 0:
                continue

            carriers_this_variant = []
            for sample_name, sample_val in zip(samples, parts[9:]):
                sample_fields = sample_val.split(":")
                if gt_idx >= len(sample_fields):
                    continue
                gt = sample_fields[gt_idx]
                z = _zygosity(gt)
                if z is None:
                    continue
                carriers_this_variant.append((sample_name, gt, z))

            for gene in plp_genes_this_variant:
                for sample_name, gt, z in carriers_this_variant:
                    stats["carrier_rows"] += 1
                    out_carriers.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{sample_name}\t{gt}\t{z}\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True, nargs="+", help="one or more fastVEP --acmg VCF outputs")
    ap.add_argument("--out-variants", required=True)
    ap.add_argument("--out-carriers", required=True)
    args = ap.parse_args()

    stats = {"total_variants": 0, "gene_rows": 0, "no_mane": 0, "multi_gene_variants": 0,
              "same_gene_multi_mane_anomaly": 0, "carrier_rows": 0}

    with open(args.out_variants, "w", encoding="utf-8") as ov, \
         open(args.out_carriers, "w", encoding="utf-8") as oc:
        ov.write(VARIANT_HEADER)
        oc.write(CARRIER_HEADER)
        for path in args.vcf:
            process_vcf(path, ov, oc, stats)

    print(f"[build_acmg_carriers] {stats['total_variants']} variants -> {stats['gene_rows']} "
          f"(variant, gene) rows ({stats['multi_gene_variants']} variants overlapped 2+ genes' "
          f"MANE transcripts); {stats['no_mane']} flagged {NO_MANE_SENTINEL}; "
          f"{stats['same_gene_multi_mane_anomaly']} same-gene multi-MANE anomalies; "
          f"{stats['carrier_rows']} carrier rows written", file=_sys.stderr)
    return 0


if __name__ == "__main__":
    _sys.exit(main())
