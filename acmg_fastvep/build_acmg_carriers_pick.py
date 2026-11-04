#!/usr/bin/env python3
"""ACMG carriers with a per-gene transcript pick, from existing fastVEP output.

Successor to build_acmg_carriers.py (MANE-only). Same inputs (already-annotated
`*.fastvep.vcf.gz`, fastVEP `--acmg`, run WITHOUT `--pick`, so every transcript's
CSQ block is present), same carrier logic, no fastVEP re-run. The difference
is which CSQ block a (variant, gene) row reports:

  * MANE-only  : the MANE_SELECT block, else the variant is flagged
                 NO_MANE_TRANSCRIPT and nothing is reported for it.
  * this script: one block per (variant, gene) chosen by the same hierarchy
                 fastVEP's `--pick-allele-gene` uses, so a gene without a MANE
                 transcript still gets a call (canonical, then TSL, ...).

The tier semantics are a line-for-line port of Huang-lab/fastVEP
crates/fastvep-annotate/src/pick.rs (`pick_score`, `appris_rank`,
`pick_key_with`) and the consequence severity table of
crates/fastvep-core/src/consequence.rs (`Consequence::rank`). Every tier
scores "lower is better"; tiers compose lexicographically in --pick-order
order; the final tie-break is the transcript ID (Feature) by plain string
order, as in fastVEP's `PickKey`.

Grouping: fastVEP groups by gene *ID* (CSQ `Gene`), not symbol, and skips
featureless "scaffold" rows (no transcript). Both are reproduced here. The
output `gene` column is still the SYMBOL, as in the MANE-only table.

Known simplification, shared with build_acmg_carriers.py: a record is treated
as one alt allele (the SpliceAI lookup and the carrier GT logic already
assume split multi-allelics). Multi-allelic ALT records are counted in the
stderr summary so this is visible if it ever happens.
"""
from __future__ import annotations

import argparse
import os as _os
import sys as _sys

_HERE = _os.path.dirname(_os.path.realpath(__file__))
_sys.path.insert(0, _HERE)
_sys.path.insert(0, _os.path.dirname(_HERE))

from plp_rules.csq import parse_csq_format, decode_csq_field  # noqa: E402
from parse_fastvep_acmg import _LABEL, _PLP, _count_criteria  # noqa: E402
from build_acmg_carriers import (  # noqa: E402
    CARRIER_HEADER,
    CSQ_FUNCTIONAL_FIELDS,
    _field,
    _is_lof,
    _open,
    _parse_spliceai_info,
    _summarize_unmapped_blocks,
    _zygosity,
)

NO_TRANSCRIPT_SENTINEL = "NO_TRANSCRIPT"

# --------------------------------------------------------------------------
# Pick hierarchy (port of fastVEP pick.rs)
# --------------------------------------------------------------------------

#: The default asked for in this pipeline: fastVEP's DEFAULT_PICK_ORDER minus
#: `appris` (the cohort's annotation source carries no APPRIS tags, so that
#: tier would score every transcript equally anyway).
DEFAULT_PICK_ORDER = "mane_select,mane_plus_clinical,canonical,tsl,biotype,ccds,rank"
#: The order used on the original Ensembl VEP path of this pipeline.
VEP_PATH_PICK_ORDER = "mane_select,canonical,tsl,biotype,rank,ccds"

_VALID_CRITERIA = ("mane_select", "mane_plus_clinical", "canonical", "appris",
                   "tsl", "biotype", "ccds", "rank")
_ALIASES = {"mane": "mane_select"}

# fastVEP `Consequence::rank` (VEP 115 ordering). Lower = more severe.
CONSEQUENCE_RANK = {
    "transcript_ablation": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "transcript_amplification": 8,
    "feature_elongation": 9,
    "feature_truncation": 10,
    "inframe_insertion": 11,
    "inframe_deletion": 12,
    "missense_variant": 13,
    "protein_altering_variant": 14,
    "splice_region_variant": 15,
    "splice_donor_5th_base_variant": 16,
    "splice_donor_region_variant": 17,
    "splice_polypyrimidine_tract_variant": 18,
    "incomplete_terminal_codon_variant": 19,
    "start_retained_variant": 20,
    "stop_retained_variant": 21,
    "synonymous_variant": 22,
    "coding_sequence_variant": 23,
    "mature_miRNA_variant": 24,
    "5_prime_UTR_variant": 25,
    "3_prime_UTR_variant": 26,
    "non_coding_transcript_exon_variant": 27,
    "intron_variant": 28,
    "NMD_transcript_variant": 29,
    "non_coding_transcript_variant": 30,
    "coding_transcript_variant": 31,
    "upstream_gene_variant": 32,
    "downstream_gene_variant": 33,
    "TFBS_ablation": 34,
    "TFBS_amplification": 35,
    "TF_binding_site_variant": 36,
    "regulatory_region_ablation": 37,
    "regulatory_region_amplification": 38,
    "regulatory_region_variant": 39,
    "intergenic_variant": 40,
    "sequence_variant": 41,
    "copy_number_change": 42,
    "copy_number_increase": 43,
    "copy_number_decrease": 44,
    "short_tandem_repeat_change": 45,
    "short_tandem_repeat_expansion": 46,
    "short_tandem_repeat_contraction": 47,
    "unidirectional_gene_fusion": 48,
    "transcript_variant": 49,
}
# fastVEP parses consequences into an enum, so it never sees an unknown term.
# Text can carry one; rank it after every known term but ahead of "no
# consequence at all" (u32::MAX in pick.rs).
_RANK_UNKNOWN_TERM = 1_000_000
_U32_MAX = 2**32 - 1
_U8_MAX = 255

FLANKING_TERMS = frozenset({"upstream_gene_variant", "downstream_gene_variant"})

# appris_rank() bands, verbatim from pick.rs.
_APPRIS_PRINCIPAL_BASE = 0
_APPRIS_PRINCIPAL_UNNUMBERED = 100
_APPRIS_ALTERNATIVE_BASE = 200
_APPRIS_ALTERNATIVE_UNNUMBERED = 300
_APPRIS_UNRECOGNISED = _U32_MAX - 1
_APPRIS_ABSENT = _U32_MAX


def parse_pick_order(spec: str) -> list[str]:
    """Mirror of pick.rs `parse_pick_order`: same names, same errors."""
    out: list[str] = []
    for raw in spec.split(","):
        name = raw.strip().lower()
        if not name:
            continue
        name = _ALIASES.get(name, name)
        if name == "length":
            raise ValueError("--pick-order: 'length' is not supported (transcript length is not "
                             "in the CSQ); remaining ties are broken by transcript ID")
        if name not in _VALID_CRITERIA:
            raise ValueError(f"--pick-order: unknown criterion {name!r}. Valid: "
                             + ", ".join(_VALID_CRITERIA))
        if name in out:
            raise ValueError(f"--pick-order: {name!r} listed more than once")
        out.append(name)
    if not out:
        raise ValueError("--pick-order: no criteria given")
    return out


def appris_rank(value: str) -> int:
    """Port of pick.rs `appris_rank`."""
    if not value:
        return _APPRIS_ABSENT
    lower = value.strip().lower()
    body = lower[len("appris_"):] if lower.startswith("appris_") else lower
    for prefix, base, unnumbered in (
        ("principal", _APPRIS_PRINCIPAL_BASE, _APPRIS_PRINCIPAL_UNNUMBERED),
        ("alternative", _APPRIS_ALTERNATIVE_BASE, _APPRIS_ALTERNATIVE_UNNUMBERED),
        ("alt", _APPRIS_ALTERNATIVE_BASE, _APPRIS_ALTERNATIVE_UNNUMBERED),
        ("p", _APPRIS_PRINCIPAL_BASE, _APPRIS_PRINCIPAL_UNNUMBERED),
        ("a", _APPRIS_ALTERNATIVE_BASE, _APPRIS_ALTERNATIVE_UNNUMBERED),
    ):
        if body.startswith(prefix):
            digits = body[len(prefix):].lstrip("_-")
            if not digits:
                return unnumbered
            if digits.isascii() and digits.isdigit():
                return base + min(max(int(digits), 1), 99)
            return unnumbered
    return _APPRIS_UNRECOGNISED


def parse_tsl(value: str) -> int:
    """fastVEP writes TSL as a bare integer (or empty). Also accept `tsl1` and
    the GFF form `1 (assigned to previous version 2)`. Absent -> u8::MAX."""
    v = value.strip().lower()
    if v.startswith("tsl"):
        v = v[3:]
    v = v.split()[0] if v.split() else ""
    if v.isdigit():
        n = int(v)
        if 0 <= n <= _U8_MAX:
            return n
    return _U8_MAX


def consequence_rank(consequence: str) -> int:
    """Min rank over '&'-joined SO terms; no terms -> u32::MAX (pick.rs Rank)."""
    ranks = [CONSEQUENCE_RANK.get(t, _RANK_UNKNOWN_TERM) for t in consequence.split("&") if t]
    return min(ranks) if ranks else _U32_MAX


def _present(v: str) -> bool:
    return bool(v) and v != "-"


class Block:
    """One CSQ block (one transcript x allele) with the fields the pick needs."""
    __slots__ = ("f", "feature", "gene_id", "symbol", "consequence",
                 "mane_select", "mane_plus_clinical", "canonical", "appris",
                 "tsl", "biotype", "ccds")

    def __init__(self, f: list[str], idx: dict[str, int]):
        g = lambda name: _field(f, idx.get(name, -1))  # noqa: E731
        self.f = f
        self.feature = g("Feature")
        self.gene_id = g("Gene")
        self.symbol = g("SYMBOL")
        self.consequence = g("Consequence")
        self.mane_select = g("MANE_SELECT")
        self.mane_plus_clinical = g("MANE_PLUS_CLINICAL")
        self.canonical = g("CANONICAL")
        self.appris = g("APPRIS")
        self.tsl = g("TSL")
        self.biotype = g("BIOTYPE")
        self.ccds = g("CCDS")

    @property
    def is_mane_select(self) -> bool:
        return _present(self.mane_select)

    @property
    def gene_key(self) -> str:
        # fastVEP groups on gene_id; fall back to the symbol for CSQ schemas
        # without a Gene subfield.
        return self.gene_id if _present(self.gene_id) else self.symbol

    def score(self, criterion: str) -> int:
        """Port of pick.rs `pick_score`. Lower is better."""
        if criterion == "mane_select":
            return 0 if _present(self.mane_select) else 1
        if criterion == "mane_plus_clinical":
            return 0 if _present(self.mane_plus_clinical) else 1
        if criterion == "canonical":
            return 0 if self.canonical.upper() == "YES" else 1
        if criterion == "appris":
            return appris_rank("" if self.appris == "-" else self.appris)
        if criterion == "tsl":
            return parse_tsl(self.tsl)
        if criterion == "biotype":
            return 0 if self.biotype == "protein_coding" else 1
        if criterion == "ccds":
            return 0 if _present(self.ccds) else 1
        if criterion == "rank":
            return consequence_rank(self.consequence)
        raise ValueError(criterion)

    def key(self, order: list[str]) -> tuple:
        return (tuple(self.score(c) for c in order), self.feature)


def pick_block(blocks: list[Block], order: list[str]) -> tuple[Block, str]:
    """Best block under `order`, plus the tier that decided it.

    The deciding tier is the first criterion at which the winner scores
    strictly better than the runner-up; `only_candidate` when the gene has a
    single block, `tie_id` when every listed criterion tied and the
    transcript-ID tie-break chose.
    """
    ranked = sorted(blocks, key=lambda b: b.key(order))
    best = ranked[0]
    if len(ranked) == 1:
        return best, "only_candidate"
    best_scores, _ = best.key(order)
    second_scores, _ = ranked[1].key(order)
    for c, a, b in zip(order, best_scores, second_scores):
        if a < b:
            return best, c
    return best, "tie_id"


def is_flanking_only(blocks: list[Block]) -> bool:
    terms = {t for b in blocks for t in b.consequence.split("&") if t}
    return bool(terms) and terms <= FLANKING_TERMS


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

VARIANT_COLUMNS = (
    "chr", "pos", "ref", "alt", "gene", "acmg_label", "acmg_criteria",
    "n_pathogenic_criteria", "n_benign_criteria", "is_acmg_PLP",
    "consequence", "impact", "hgvsc", "hgvsp", "exon", "intron", "is_lof_consequence",
    "spliceai_ds_ag", "spliceai_ds_al", "spliceai_ds_dg", "spliceai_ds_dl",
    "transcript_id", "pick_tier", "is_mane_select", "acmg_label_mane", "mane_discordant",
)
VARIANT_HEADER = "\t".join(VARIANT_COLUMNS) + "\n"

_REQUIRED_FIELDS = ("ACMG", "Feature", "Consequence")
_PICK_FIELDS = {
    "mane_select": "MANE_SELECT", "mane_plus_clinical": "MANE_PLUS_CLINICAL",
    "canonical": "CANONICAL", "appris": "APPRIS", "tsl": "TSL",
    "biotype": "BIOTYPE", "ccds": "CCDS", "rank": "Consequence",
}


def new_stats() -> dict:
    return {k: 0 for k in (
        "total_variants", "gene_rows", "multi_gene_variants", "no_transcript",
        "flanking_rows_dropped", "flanking_plp_rows_dropped", "variants_all_flanking",
        "mane_rows", "fallback_rows", "mane_discordant_rows", "plp_rows",
        "multiallelic_records", "carrier_rows", "hemi_calls",
    )} | {"by_tier": {}}


def _acmg_of(f, idx):
    sh = _field(f, idx["ACMG"])
    crit = _field(f, idx.get("ACMG_CRITERIA", -1))
    codes = [c for c in crit.replace("&", ";").split(";") if c]
    return sh, codes


def process_vcf(path: str, out_variants, out_carriers, stats: dict, order: list[str],
                drop_flanking: bool = True) -> None:
    idx: dict[str, int] = {}
    spliceai_id = None
    samples: list[str] = []
    gt_idx_cache: dict[str, int] = {}
    have_schema = False

    with _open(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=SpliceAI,"):
                spliceai_id = "SpliceAI"
                continue
            if line.startswith("##INFO=<ID=CSQ,"):
                schema = parse_csq_format(line)
                idx = {name: i for i, name in enumerate(schema.fields)}
                missing = [f for f in _REQUIRED_FIELDS if f not in idx]
                if missing:
                    print(f"ERROR: CSQ ({path}) lacks {missing} — was fastVEP run with --acmg?",
                          file=_sys.stderr)
                    _sys.exit(2)
                absent = sorted({_PICK_FIELDS[c] for c in order} - set(idx))
                if absent:
                    # A tier whose field is absent scores every transcript the
                    # same, so it silently decides nothing. Say so once per file.
                    print(f"WARN: CSQ ({path}) has no {absent} — those --pick-order tiers "
                          "cannot separate transcripts", file=_sys.stderr)
                have_schema = True
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:] if len(cols) > 9 else []
                continue
            if line.startswith("#"):
                continue
            if not have_schema:
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
            if "," in alt:
                stats["multiallelic_records"] += 1
            raw_blocks = csq.split(",")

            # Group candidate blocks by gene. Scaffold / featureless blocks
            # (no transcript) and blocks with no gene (intergenic, regulatory)
            # are never candidates, as in pick.rs.
            by_gene: dict[str, list[Block]] = {}
            for raw in raw_blocks:
                b = Block(raw.split("|"), idx)
                if not _present(b.feature) or not _present(b.gene_key):
                    continue
                by_gene.setdefault(b.gene_key, []).append(b)

            if not by_gene:
                stats["no_transcript"] += 1
                genes_seen, cons_seen = _summarize_unmapped_blocks(
                    raw_blocks, idx.get("SYMBOL", -1), idx["Consequence"])
                row = [chrom, pos, ref, alt, genes_seen, NO_TRANSCRIPT_SENTINEL, "", "0", "0", "0",
                       cons_seen] + [""] * (len(VARIANT_COLUMNS) - 11)
                out_variants.write("\t".join(row) + "\n")
                continue

            spliceai_by_allele = _parse_spliceai_info(info, spliceai_id)
            ds = spliceai_by_allele.get(alt, ("", "", "", ""))

            rows_written = 0
            plp_genes: list[str] = []
            for gene_blocks in by_gene.values():
                picked, tier = pick_block(gene_blocks, order)
                if drop_flanking and is_flanking_only(gene_blocks):
                    stats["flanking_rows_dropped"] += 1
                    if _field(picked.f, idx["ACMG"]) in _PLP:
                        stats["flanking_plp_rows_dropped"] += 1
                    continue

                sh, codes = _acmg_of(picked.f, idx)
                npath, nben = _count_criteria(codes)
                is_plp = int(sh in _PLP)
                func = {n: _field(picked.f, idx.get(n, -1)) for n in CSQ_FUNCTIONAL_FIELDS}

                mane_blocks = sorted((b for b in gene_blocks if b.is_mane_select),
                                     key=lambda b: b.feature)
                if mane_blocks:
                    m_sh, _ = _acmg_of(mane_blocks[0].f, idx)
                    label_mane = _LABEL.get(m_sh, m_sh)
                    discordant = int(is_plp != int(m_sh in _PLP))
                    stats["mane_discordant_rows"] += discordant
                    discordant_s = str(discordant)
                else:
                    label_mane, discordant_s = "", ""

                gene = picked.symbol or picked.gene_id
                row = [chrom, pos, ref, alt, gene, _LABEL.get(sh, sh), ";".join(codes),
                       str(npath), str(nben), str(is_plp),
                       func["Consequence"], func["IMPACT"], func["HGVSc"], func["HGVSp"],
                       func["EXON"], func["INTRON"], str(_is_lof(func["Consequence"])),
                       ds[0], ds[1], ds[2], ds[3],
                       picked.feature, tier, str(int(picked.is_mane_select)),
                       label_mane, discordant_s]
                out_variants.write("\t".join(row) + "\n")
                rows_written += 1
                stats["gene_rows"] += 1
                stats["by_tier"][tier] = stats["by_tier"].get(tier, 0) + 1
                stats["mane_rows" if picked.is_mane_select else "fallback_rows"] += 1
                if is_plp:
                    stats["plp_rows"] += 1
                    plp_genes.append(gene)

            if rows_written > 1:
                stats["multi_gene_variants"] += 1
            if rows_written == 0:
                stats["variants_all_flanking"] += 1

            if not plp_genes or not samples:
                continue

            format_field = parts[8] if len(parts) > 8 else ""
            if format_field not in gt_idx_cache:
                gt_fields = format_field.split(":")
                gt_idx_cache[format_field] = gt_fields.index("GT") if "GT" in gt_fields else -1
            gt_idx = gt_idx_cache[format_field]
            if gt_idx < 0:
                continue

            carriers = []
            for sample_name, sample_val in zip(samples, parts[9:]):
                sf = sample_val.split(":")
                if gt_idx >= len(sf):
                    continue
                z = _zygosity(sf[gt_idx])
                if z is None:
                    continue
                if z == "hemi":
                    stats["hemi_calls"] += 1
                carriers.append((sample_name, sf[gt_idx], z))

            for gene in plp_genes:
                for sample_name, gt, z in carriers:
                    stats["carrier_rows"] += 1
                    out_carriers.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{sample_name}\t{gt}\t{z}\n")


def format_summary(stats: dict, order: list[str], drop_flanking: bool) -> str:
    tiers = ", ".join(f"{k}={v}" for k, v in sorted(stats["by_tier"].items(),
                                                       key=lambda kv: -kv[1]))
    return (
        f"[build_acmg_carriers_pick] order={','.join(order)} drop_flanking={int(drop_flanking)}\n"
        f"  {stats['total_variants']} variants -> {stats['gene_rows']} (variant, gene) rows "
        f"({stats['multi_gene_variants']} variants with 2+ gene rows)\n"
        f"  picked MANE Select: {stats['mane_rows']}; fallback (non-MANE): {stats['fallback_rows']}\n"
        f"  pick_tier: {tiers or '-'}\n"
        f"  flanking-only gene rows dropped: {stats['flanking_rows_dropped']} "
        f"(of which P/LP: {stats['flanking_plp_rows_dropped']}); "
        f"variants left with no row because every gene was flanking-only: "
        f"{stats['variants_all_flanking']}\n"
        f"  {NO_TRANSCRIPT_SENTINEL} (no gene at all): {stats['no_transcript']}\n"
        f"  P/LP rows: {stats['plp_rows']}; mane_discordant rows: {stats['mane_discordant_rows']}\n"
        f"  multi-allelic ALT records: {stats['multiallelic_records']}\n"
        f"  carrier rows: {stats['carrier_rows']} (hemizygous calls: {stats['hemi_calls']})"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, nargs="+", help="fastVEP --acmg VCF output(s), no --pick")
    ap.add_argument("--out-variants", required=True)
    ap.add_argument("--out-carriers", required=True)
    ap.add_argument("--pick-order", default=DEFAULT_PICK_ORDER,
                    help=f"comma-separated tiers (default {DEFAULT_PICK_ORDER}; VEP-path order: "
                         f"{VEP_PATH_PICK_ORDER}). Final tie-break is always transcript ID.")
    fl = ap.add_mutually_exclusive_group()
    fl.add_argument("--drop-flanking-genes", dest="drop_flanking", action="store_true",
                    default=True, help="drop gene rows whose only consequences are "
                    "upstream/downstream_gene_variant (default)")
    fl.add_argument("--keep-flanking-genes", dest="drop_flanking", action="store_false")
    args = ap.parse_args(argv)

    try:
        order = parse_pick_order(args.pick_order)
    except ValueError as e:
        ap.error(str(e))

    stats = new_stats()
    with open(args.out_variants, "w", encoding="utf-8") as ov, \
         open(args.out_carriers, "w", encoding="utf-8") as oc:
        ov.write(VARIANT_HEADER)
        oc.write(CARRIER_HEADER)
        for path in args.vcf:
            process_vcf(path, ov, oc, stats, order, args.drop_flanking)

    print(format_summary(stats, order, args.drop_flanking), file=_sys.stderr)
    return 0


if __name__ == "__main__":
    _sys.exit(main())
