#!/usr/bin/env python3
"""ACMG carriers with a per-gene transcript pick, from existing fastVEP output.

Successor to build_acmg_carriers.py (MANE-only). Same inputs (already-annotated
`*.fastvep.vcf.gz`, fastVEP `--acmg`, run WITHOUT `--pick`, so every transcript's
CSQ block is present), same carrier logic, no fastVEP re-run. The difference
is which CSQ block a (variant, gene) row reports:

  * MANE-only  : the MANE_SELECT block, else the variant is flagged
                 NO_MANE_TRANSCRIPT and nothing is reported for it.
  * this script: one block per (variant, gene), so a gene without a MANE
                 transcript still gets a call (canonical, then TSL, ...).

Two ways of choosing that block (--selection):

  * fastvep_pick (DEFAULT; fastVEP v0.5.0 output from run_fastvep.sh --acmg,
    which passes --flag-pick-allele-gene): the block fastVEP itself marked
    PICK=1 for that (allele, gene). No pick logic is re-implemented here; the
    CSQ must carry a PICK field (error otherwise). An (allele, gene) with no
    PICK=1 block, or with more than one, is not resolved: its row is flagged
    (acmg_label NO_PICK_FLAG / MULTIPLE_PICK_FLAGS, is_acmg_PLP 0, pick_flag
    says which transcripts) and counted in the summary. acmg_label_mane is the
    MANE Select block's label, or NO_MANE_TRANSCRIPT. Per fastVEP's pick.rs,
    --flag-pick-allele-gene flags exactly one block per (allele, gene), so a
    flagged row means this parser and fastVEP's output disagree: both tables
    are still written, then the script exits 3 (--allow-pick-anomalies accepts).
  * python_pick (DEPRECATED; kept for a one-time comparison against the v6
    outputs, which have no PICK field): the Python port of fastVEP's pick
    hierarchy described below.

Both read fields by CSQ header name, so the trailing ACMG_DISORDERS field of
v0.5.0 (`INHERITANCE:CLASS:CRITERIA` per mode of inheritance, joined by '&')
is tolerated and copied to the acmg_disorders column.

python_pick: the tier semantics are a line-for-line port of Huang-lab/fastVEP
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
import json
import os as _os
import sys as _sys

_HERE = _os.path.dirname(_os.path.realpath(__file__))
_sys.path.insert(0, _HERE)
_sys.path.insert(0, _os.path.dirname(_HERE))

from plp_rules.csq import parse_csq_format, decode_csq_field  # noqa: E402
import carrier_markers  # noqa: E402
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
NO_MANE_SENTINEL = "NO_MANE_TRANSCRIPT"
NO_PICK_SENTINEL = "NO_PICK_FLAG"
MULTI_PICK_SENTINEL = "MULTIPLE_PICK_FLAGS"
SELECTIONS = ("fastvep_pick", "python_pick")

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
                 "tsl", "biotype", "ccds", "allele", "pick")

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
        self.allele = g("Allele")
        self.pick = g("PICK")      # "1" on fastVEP's own choice (--flag-pick-allele-gene)

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
    # picked transcript (CSQ)
    "protein_position", "amino_acids", "cds_position", "flags", "existing_variation",
    "is_last_exon",
    # variant / gene level (fastVEP FV_* INFO fields)
    "clinvar_significance", "clinvar_review_status",
    "gnomad_af", "gnomad_grpmax_af", "gnomad_grpmax_pop", "gnomad_hom_count",
    "revel", "gene_pli", "gene_loeuf",
    # fastVEP v0.5.0: per-inheritance-mode calls (ACMG_DISORDERS), and whether
    # fastVEP's PICK flag resolved this (allele, gene) to exactly one block
    "acmg_disorders", "pick_flag",
)
VARIANT_HEADER = "\t".join(VARIANT_COLUMNS) + "\n"

_REQUIRED_FIELDS = ("ACMG", "Feature", "Consequence")
_PICK_FIELDS = {
    "mane_select": "MANE_SELECT", "mane_plus_clinical": "MANE_PLUS_CLINICAL",
    "canonical": "CANONICAL", "appris": "APPRIS", "tsl": "TSL",
    "biotype": "BIOTYPE", "ccds": "CCDS", "rank": "Consequence",
}


# Extra per-transcript CSQ fields copied from the picked block.
CSQ_EXTRA_FIELDS = ("Protein_position", "Amino_acids", "CDS_position", "FLAGS",
                    "Existing_variation")
# Allele-keyed / gene-keyed INFO fields fastVEP writes from its SA databases.
FV_INFO_IDS = ("FV_CLINVAR", "FV_GNOMAD", "FV_REVEL", "FV_GNOMAD_GENE")
# gnomAD v4 "grpmax": highest AF among non-bottlenecked genetic ancestry
# groups, i.e. excluding Amish, Ashkenazi Jewish, Finnish, Middle Eastern and
# remaining/other. This is the frequency the ACMG BA1/BS1 rules refer to.
GNOMAD_GRPMAX_POPS = ("AFR", "AMR", "EAS", "NFE", "SAS")


def vep_allele(ref: str, alt: str) -> str:
    """The VEP/fastVEP spelling of an alt allele: for an indel whose REF and ALT
    share a first base, that anchor base is dropped ('-' if nothing is left)."""
    if len(ref) != len(alt) and ref and alt and ref[0] == alt[0]:
        return alt[1:] or "-"
    return alt


def is_last_exon(exon: str) -> str:
    """'1' if the CSQ EXON value ('N/M' or 'N-K/M') reaches the last exon, '0'
    if not, '' when there is no exon number (intronic, splice-site only...)."""
    if "/" not in exon:
        return ""
    num, _, total = exon.partition("/")
    try:
        last = max(int(x) for x in num.split("-") if x)
        return str(int(last == int(total)))
    except ValueError:
        return ""


def _fv_entries(info_kv: dict, fid: str, formats: dict) -> list[dict]:
    val = info_kv.get(fid)
    if not val or fid not in formats:
        return []
    fields = formats[fid]
    return [dict(zip(fields, (decode_csq_field(x) for x in e.split("|"))))
            for e in val.split(",")]


def _float(v: str):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def fv_variant_columns(info_kv: dict, formats: dict, allele_keys: set, stats: dict):
    """(clinvar_sig, clinvar_review, gnomad_af, grpmax_af, grpmax_pop, hom_count,
    revel) for this record's alt allele. Empty strings where absent."""
    out = {}
    unmatched = False
    for fid in ("FV_CLINVAR", "FV_GNOMAD", "FV_REVEL"):
        entries = _fv_entries(info_kv, fid, formats)
        hit = [e for e in entries if e.get("ALLELE", "") in allele_keys]
        if entries and not hit:
            unmatched = True
        out[fid] = hit
    if unmatched:
        stats["fv_allele_unmatched"] += 1

    def uniq(entries, key):
        seen = []
        for e in entries:
            v = e.get(key, "")
            if v and v not in seen:
                seen.append(v)
        return ";".join(seen)

    sig = uniq(out["FV_CLINVAR"], "SIGNIFICANCE")
    review = uniq(out["FV_CLINVAR"], "REVIEW_STATUS")

    af = grp_af = grp_pop = hom = ""
    if out["FV_GNOMAD"]:
        g = out["FV_GNOMAD"][0]
        af, hom = g.get("ALL_AF", ""), g.get("ALL_HC", "")
        best = None
        for pop in GNOMAD_GRPMAX_POPS:
            x = _float(g.get(f"{pop}_AF", ""))
            if x is not None and (best is None or x > best[0]):
                best = (x, pop, g[f"{pop}_AF"])
        if best is not None:
            grp_af, grp_pop = best[2], best[1]

    revels = [x for x in (_float(e.get("SCORE")) for e in out["FV_REVEL"]) if x is not None]
    revel = f"{max(revels):g}" if revels else ""
    return [sig, review, af, grp_af, grp_pop, hom, revel]


def fv_gene_columns(info_kv: dict, formats: dict, symbol: str):
    """(pLI, LOEUF) for this gene symbol from FV_GNOMAD_GENE."""
    for e in _fv_entries(info_kv, "FV_GNOMAD_GENE", formats):
        if e.get("SYMBOL") == symbol:
            return [e.get("PLI", ""), e.get("LOEUF", "")]
    return ["", ""]


def new_stats() -> dict:
    return {k: 0 for k in (
        "total_variants", "gene_rows", "multi_gene_variants", "no_transcript",
        "flanking_rows_dropped", "flanking_plp_rows_dropped", "variants_all_flanking",
        "mane_rows", "fallback_rows", "mane_discordant_rows", "plp_rows",
        "multiallelic_records", "carrier_rows", "hemi_calls", "fv_allele_unmatched",
        "no_pick_rows", "multi_pick_rows", "no_mane_rows",
    )} | {"by_tier": {}}


def fastvep_pick(blocks: list[Block], idx: dict[str, int]):
    """fastVEP's own choice for one (allele, gene): the block flagged PICK=1.

    Returns (block, tier, pick_flag). block is None when fastVEP flagged none
    or several; nothing is chosen here in that case."""
    picks = [b for b in blocks if b.pick == "1"]
    if len(picks) == 1:
        return picks[0], "fastvep_pick", "ok"
    if not picks:
        return None, "", "no_pick"
    detail = ",".join(f"{b.feature}={_LABEL.get(_acmg_of(b.f, idx)[0], _acmg_of(b.f, idx)[0])}"
                      for b in sorted(picks, key=lambda b: b.feature))
    return None, "", f"multi_pick:{detail}"


_MAX_PICK_WARNINGS = 20


def _acmg_of(f, idx):
    sh = _field(f, idx["ACMG"])
    crit = _field(f, idx.get("ACMG_CRITERIA", -1))
    codes = [c for c in crit.replace("&", ";").split(";") if c]
    return sh, codes


def process_vcf(path: str, out_variants, out_carriers, stats: dict, order: list[str],
                drop_flanking: bool = True, selection: str = "fastvep_pick") -> None:
    idx: dict[str, int] = {}
    spliceai_id = None
    fv_formats: dict[str, tuple] = {}
    samples: list[str] = []
    gt_idx_cache: dict[str, int] = {}
    have_schema = False

    with _open(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=SpliceAI,"):
                spliceai_id = "SpliceAI"
                continue
            if line.startswith("##INFO=<ID=FV_"):
                fid = line[len("##INFO=<ID="):].split(",", 1)[0]
                if fid in FV_INFO_IDS:
                    try:
                        fv_formats[fid] = parse_csq_format(line).fields
                    except ValueError:
                        pass
                continue
            if line.startswith("##INFO=<ID=CSQ,"):
                schema = parse_csq_format(line)
                idx = {name: i for i, name in enumerate(schema.fields)}
                missing = [f for f in _REQUIRED_FIELDS if f not in idx]
                if missing:
                    print(f"ERROR: CSQ ({path}) lacks {missing} — was fastVEP run with --acmg?",
                          file=_sys.stderr)
                    _sys.exit(2)
                if selection == "fastvep_pick" and "PICK" not in idx:
                    print(f"ERROR: CSQ ({path}) has no PICK field. --selection fastvep_pick reads "
                          "fastVEP's own PICK=1 flag, written only by fastVEP >= v0.5.0 run with "
                          "--flag-pick-allele-gene (run_fastvep.sh --acmg). For older output use "
                          "--selection python_pick (deprecated) or build_acmg_carriers.py (MANE).",
                          file=_sys.stderr)
                    _sys.exit(2)
                absent = sorted({_PICK_FIELDS[c] for c in order} - set(idx))
                if absent and selection == "python_pick":
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
            info_kv = dict(kv.split("=", 1) if "=" in kv else (kv, "")
                           for kv in info.split(";"))
            csq = info_kv.get("CSQ", "")
            if not csq:
                continue

            stats["total_variants"] += 1
            if "," in alt:
                stats["multiallelic_records"] += 1
            raw_blocks = csq.split(",")

            # Group candidate blocks by gene. Scaffold / featureless blocks
            # (no transcript) and blocks with no gene (intergenic, regulatory)
            # are never candidates, as in pick.rs.
            # fastvep_pick groups per (allele, gene), as --flag-pick-allele-gene does.
            by_gene: dict[tuple, list[Block]] = {}
            for raw in raw_blocks:
                b = Block(raw.split("|"), idx)
                if not _present(b.feature) or not _present(b.gene_key):
                    continue
                key = (b.allele, b.gene_key) if selection == "fastvep_pick" else ("", b.gene_key)
                by_gene.setdefault(key, []).append(b)

            if not by_gene:
                stats["no_transcript"] += 1
                genes_seen, cons_seen = _summarize_unmapped_blocks(
                    raw_blocks, idx.get("SYMBOL", -1), idx["Consequence"])
                row = [chrom, pos, ref, alt, genes_seen, NO_TRANSCRIPT_SENTINEL, "", "0", "0", "0",
                       cons_seen] + [""] * (len(VARIANT_COLUMNS) - 11)
                out_variants.write("\t".join(row) + "\n")
                continue

            allele_keys = {alt, vep_allele(ref, alt)}
            spliceai_by_allele = _parse_spliceai_info(info, spliceai_id)
            ds = next((spliceai_by_allele[a] for a in (alt, vep_allele(ref, alt))
                       if a in spliceai_by_allele), ("", "", "", ""))
            fv_variant = fv_variant_columns(info_kv, fv_formats, allele_keys, stats)

            rows_written = 0
            plp_genes: list[str] = []
            for gene_blocks in by_gene.values():
                if selection == "fastvep_pick":
                    picked, tier, pick_flag = fastvep_pick(gene_blocks, idx)
                else:
                    (picked, tier), pick_flag = pick_block(gene_blocks, order), ""
                if drop_flanking and is_flanking_only(gene_blocks):
                    stats["flanking_rows_dropped"] += 1
                    if picked is not None and _field(picked.f, idx["ACMG"]) in _PLP:
                        stats["flanking_plp_rows_dropped"] += 1
                    continue

                mane_blocks = sorted((b for b in gene_blocks if b.is_mane_select),
                                     key=lambda b: b.feature)
                m_sh = _acmg_of(mane_blocks[0].f, idx)[0] if mane_blocks else None

                if picked is None:
                    # fastVEP flagged no block, or several, for this (allele, gene).
                    # Reported, never resolved here: that would be a pick of our own.
                    first = gene_blocks[0]
                    sentinel = NO_PICK_SENTINEL if pick_flag == "no_pick" else MULTI_PICK_SENTINEL
                    stats["no_pick_rows" if pick_flag == "no_pick" else "multi_pick_rows"] += 1
                    n_warn = stats["no_pick_rows"] + stats["multi_pick_rows"]
                    if n_warn <= _MAX_PICK_WARNINGS:
                        print(f"WARN: {chrom}:{pos} {ref}>{alt} gene={first.symbol or first.gene_id}: "
                              f"{pick_flag} (fastVEP PICK=1 on {len([b for b in gene_blocks if b.pick == '1'])} "
                              f"of {len(gene_blocks)} blocks); row flagged {sentinel}, not P/LP",
                              file=_sys.stderr)
                    picks = [b.feature for b in gene_blocks if b.pick == "1"]
                    row = [chrom, pos, ref, alt, first.symbol or first.gene_id, sentinel, "",
                           "0", "0", "0",
                           ";".join(sorted({b.consequence for b in gene_blocks})),
                           "", "", "", "", "", "", "", "", "", "",
                           ";".join(picks), "", "",
                           _LABEL.get(m_sh, m_sh) if m_sh is not None else NO_MANE_SENTINEL, "",
                           "", "", "", "", "", "",
                           *fv_variant, *fv_gene_columns(info_kv, fv_formats, first.symbol),
                           "", pick_flag]
                    out_variants.write("\t".join(row) + "\n")
                    rows_written += 1
                    stats["gene_rows"] += 1
                    continue

                sh, codes = _acmg_of(picked.f, idx)
                npath, nben = _count_criteria(codes)
                is_plp = int(sh in _PLP)
                func = {n: _field(picked.f, idx.get(n, -1)) for n in CSQ_FUNCTIONAL_FIELDS}

                if m_sh is not None:
                    label_mane = _LABEL.get(m_sh, m_sh)
                    discordant = int(is_plp != int(m_sh in _PLP))
                    stats["mane_discordant_rows"] += discordant
                    discordant_s = str(discordant)
                else:
                    # python_pick keeps its v6 output (blank); fastvep_pick says so.
                    label_mane = NO_MANE_SENTINEL if selection == "fastvep_pick" else ""
                    discordant_s = ""
                    stats["no_mane_rows"] += 1

                gene = picked.symbol or picked.gene_id
                row = [chrom, pos, ref, alt, gene, _LABEL.get(sh, sh), ";".join(codes),
                       str(npath), str(nben), str(is_plp),
                       func["Consequence"], func["IMPACT"], func["HGVSc"], func["HGVSp"],
                       func["EXON"], func["INTRON"], str(_is_lof(func["Consequence"])),
                       ds[0], ds[1], ds[2], ds[3],
                       picked.feature, tier, str(int(picked.is_mane_select)),
                       label_mane, discordant_s,
                       *(_field(picked.f, idx.get(n, -1)) for n in CSQ_EXTRA_FIELDS),
                       is_last_exon(func["EXON"]),
                       *fv_variant, *fv_gene_columns(info_kv, fv_formats, picked.symbol),
                       _field(picked.f, idx.get("ACMG_DISORDERS", -1)), pick_flag]
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


def format_summary(stats: dict, order: list[str], drop_flanking: bool,
                   selection: str = "python_pick") -> str:
    tiers = ", ".join(f"{k}={v}" for k, v in sorted(stats["by_tier"].items(),
                                                       key=lambda kv: -kv[1]))
    how = ("selection=fastvep_pick (fastVEP PICK=1)" if selection == "fastvep_pick"
           else f"selection=python_pick (DEPRECATED) order={','.join(order)}")
    return (
        f"[build_acmg_carriers_pick] {how} drop_flanking={int(drop_flanking)}\n"
        f"  fastVEP PICK flag: (allele, gene) rows with no PICK=1: {stats.get('no_pick_rows', 0)}; "
        f"with more than one: {stats.get('multi_pick_rows', 0)}"
        f"{'  <-- FLAGGED ROWS, see stderr' if stats.get('no_pick_rows', 0) + stats.get('multi_pick_rows', 0) else ''}\n"
        f"  rows whose gene has no MANE Select block: {stats.get('no_mane_rows', 0)}\n"
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
        f"  multi-allelic ALT records: {stats['multiallelic_records']}; variants whose "
        f"FV_CLINVAR/FV_GNOMAD/FV_REVEL had no entry for the alt allele: "
        f"{stats['fv_allele_unmatched']}\n"
        f"  carrier rows: {stats['carrier_rows']} (hemizygous calls: {stats['hemi_calls']})"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, nargs="+", help="fastVEP --acmg VCF output(s), no --pick")
    ap.add_argument("--out-variants", required=True)
    ap.add_argument("--out-carriers", required=True)
    ap.add_argument("--selection", choices=SELECTIONS, default="fastvep_pick",
                    help="fastvep_pick (default): fastVEP's own PICK=1 block (v0.5.0 output, "
                         "--flag-pick-allele-gene). python_pick: DEPRECATED Python port of the "
                         "pick order, for comparing against PICK-less v6 output only.")
    ap.add_argument("--pick-order", default=None,
                    help=f"python_pick only: comma-separated tiers (default {DEFAULT_PICK_ORDER}; "
                         f"VEP-path order: {VEP_PATH_PICK_ORDER}). Final tie-break is always "
                         "transcript ID.")
    fl = ap.add_mutually_exclusive_group()
    fl.add_argument("--drop-flanking-genes", dest="drop_flanking", action="store_true",
                    default=True, help="drop gene rows whose only consequences are "
                    "upstream/downstream_gene_variant (default)")
    fl.add_argument("--keep-flanking-genes", dest="drop_flanking", action="store_false")
    ap.add_argument("--stats-json", help="also write the summary counters here (for per-chromosome "
                    "aggregation by run_acmg_carriers_chr.sh)")
    ap.add_argument("--no-markers", action="store_true",
                    help="do not write DONE/FAILED next to --out-variants (run_acmg_carriers_chr.sh "
                         "writes them for the whole chromosome instead)")
    ap.add_argument("--allow-pick-anomalies", action="store_true",
                    help="fastvep_pick: exit 0 even when an (allele, gene) had no or several "
                         "PICK=1 blocks (default: exit 3 after writing the outputs)")
    args = ap.parse_args(argv)

    if args.selection == "fastvep_pick" and args.pick_order is not None:
        ap.error("--pick-order applies to --selection python_pick only; fastvep_pick uses "
                 "fastVEP's own PICK flag")
    try:
        order = parse_pick_order(args.pick_order or DEFAULT_PICK_ORDER)
    except ValueError as e:
        ap.error(str(e))
    if args.selection == "python_pick":
        print("WARN: --selection python_pick is DEPRECATED: it re-implements fastVEP's pick "
              "order in Python. Use it only to compare against PICK-less (v6) output; fastVEP "
              "v0.5.0 output should use the default fastvep_pick.", file=_sys.stderr)

    stats = new_stats()
    with open(args.out_variants, "w", encoding="utf-8") as ov, \
         open(args.out_carriers, "w", encoding="utf-8") as oc:
        ov.write(VARIANT_HEADER)
        oc.write(CARRIER_HEADER)
        for path in args.vcf:
            process_vcf(path, ov, oc, stats, order, args.drop_flanking, args.selection)

    print(format_summary(stats, order, args.drop_flanking, args.selection), file=_sys.stderr)
    if args.stats_json:
        with open(args.stats_json, "w", encoding="utf-8") as fh:
            json.dump({"order": order, "drop_flanking": args.drop_flanking,
                       "selection": args.selection, "stats": stats}, fh)
    rc = _anomaly_exit(stats, args.selection, args.allow_pick_anomalies)
    if not args.no_markers:
        d = _os.path.dirname(_os.path.abspath(args.out_variants))
        counts = (f"no PICK=1: {stats.get('no_pick_rows', 0)}; "
                  f"more than one: {stats.get('multi_pick_rows', 0)}")
        if rc:
            carrier_markers.write_failed(d, "PICK anomalies: (allele, gene) rows without exactly "
                                            "one PICK=1", rc, counts)
        else:
            carrier_markers.write_done(d, f"selection={args.selection} "
                                          f"allow_pick_anomalies={int(args.allow_pick_anomalies)} {counts}")
    return rc


def _anomaly_exit(stats: dict, selection: str, allow: bool) -> int:
    n = stats.get("no_pick_rows", 0) + stats.get("multi_pick_rows", 0)
    if selection != "fastvep_pick" or n == 0:
        return 0
    if allow:
        print(f"WARN: {n} (allele, gene) row(s) without exactly one PICK=1 block; accepted "
              "(--allow-pick-anomalies)", file=_sys.stderr)
        return 0
    print(f"ERROR: {n} (allele, gene) row(s) had no or several PICK=1 blocks. fastVEP's "
          "--flag-pick-allele-gene flags exactly one, so this output and this parser disagree "
          "(outputs were written; rows flagged NO_PICK_FLAG / MULTIPLE_PICK_FLAGS). Investigate, "
          "or pass --allow-pick-anomalies.", file=_sys.stderr)
    return 3


def aggregate_stats_main(argv=None) -> int:
    """`build_acmg_carriers_pick.py --aggregate a.json b.json ...`: sum per-chunk
    counters and print the same summary for the whole set."""
    paths = list(argv if argv is not None else _sys.argv[2:])
    allow = "--allow-pick-anomalies" in paths
    paths = [p for p in paths if p != "--allow-pick-anomalies"]
    total = new_stats()
    order, drop, selection = [], True, "python_pick"
    for p in paths:
        with open(p, encoding="utf-8") as fh:
            d = json.load(fh)
        order, drop = d["order"], d["drop_flanking"]
        selection = d.get("selection", "python_pick")
        for k, v in d["stats"].items():
            if k == "by_tier":
                for t, n in v.items():
                    total["by_tier"][t] = total["by_tier"].get(t, 0) + n
            else:
                total[k] = total.get(k, 0) + v
    print(format_summary(total, order, drop, selection))
    return _anomaly_exit(total, selection, allow)


if __name__ == "__main__":
    if len(_sys.argv) > 1 and _sys.argv[1] == "--aggregate":
        _sys.exit(aggregate_stats_main())
    _sys.exit(main())
