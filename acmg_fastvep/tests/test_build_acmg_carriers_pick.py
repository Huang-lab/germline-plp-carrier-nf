"""Synthetic-CSQ tests for build_acmg_carriers_pick.py (per-gene pick fallback)."""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent.parent))

import build_acmg_carriers_pick as P  # noqa: E402
from build_acmg_carriers import _zygosity  # noqa: E402

SCRIPT = HERE.parent / "build_acmg_carriers_pick.py"

FIELDS = ("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|"
          "CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|CCDS|ACMG|ACMG_CRITERIA")
CSQ_HEADER = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from '
              f'fastVEP. Format: {FIELDS}">')


def blk(cons, sym, gene, feat, acmg="", crit="", *, biotype="protein_coding", canonical="",
        mane="", mpc="", tsl="", ccds="", impact="MODERATE", allele="T"):
    return "|".join([allele, cons, impact, sym, gene, feat, biotype, "", "", "", "",
                     canonical, mane, mpc, tsl, ccds, acmg, crit])


def write_vcf(tmp_path, records, samples=("S1", "S2", "S3")):
    lines = ["##fileformat=VCFv4.2", CSQ_HEADER,
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                        *samples])]
    for chrom, pos, blocks, gts in records:
        lines.append("\t".join([chrom, str(pos), ".", "C", "T", ".", ".",
                                "CSQ=" + ",".join(blocks), "GT", *gts]))
    p = tmp_path / "in.fastvep.vcf"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def run(tmp_path, records, *extra):
    vcf = write_vcf(tmp_path, records)
    ov, oc = tmp_path / "v.tsv", tmp_path / "c.tsv"
    res = subprocess.run([sys.executable, str(SCRIPT), "--vcf", str(vcf), "--out-variants", str(ov),
                          "--out-carriers", str(oc), *extra],
                         check=True, capture_output=True, text=True)
    variants = list(csv.DictReader(ov.open(), delimiter="\t"))
    carriers = list(csv.DictReader(oc.open(), delimiter="\t"))
    return variants, carriers, res.stderr


HOMREF = ("0/0", "0/0", "0/0")


def test_mane_present_is_picked_over_more_severe_non_mane(tmp_path):
    rec = ("17", 100, [
        blk("missense_variant", "BRCA1", "G1", "ENST2", "VUS", "PM2", canonical="YES",
            mane="NM_1.1", tsl="1"),
        blk("stop_gained", "BRCA1", "G1", "ENST1", "P", "PVS1&PM2", tsl="1", impact="HIGH"),
    ], ("0/1", "0/0", "0/0"))
    v, c, _ = run(tmp_path, [rec])
    assert len(v) == 1
    r = v[0]
    assert (r["transcript_id"], r["pick_tier"], r["is_mane_select"]) == ("ENST2", "mane_select", "1")
    assert r["acmg_label"] == "Uncertain_significance"
    assert r["acmg_label_mane"] == "Uncertain_significance" and r["mane_discordant"] == "0"
    assert c == []

    # With rank first, the non-MANE stop_gained wins and is flagged discordant.
    v, c, _ = run(tmp_path, [rec], "--pick-order", "rank,mane_select")
    r = v[0]
    assert (r["transcript_id"], r["pick_tier"], r["is_mane_select"]) == ("ENST1", "rank", "0")
    assert r["acmg_label"] == "Pathogenic" and r["mane_discordant"] == "1"
    assert r["acmg_label_mane"] == "Uncertain_significance"
    assert [x["person_id"] for x in c] == ["S1"]


def test_no_mane_falls_back_to_canonical(tmp_path):
    rec = ("1", 200, [
        blk("missense_variant", "GENEX", "G2", "ENST10", "VUS", "", tsl="1"),
        blk("stop_gained", "GENEX", "G2", "ENST11", "LP", "PVS1", canonical="YES", tsl="2"),
    ], HOMREF)
    v, _, _ = run(tmp_path, [rec])
    r = v[0]
    assert (r["transcript_id"], r["pick_tier"], r["is_mane_select"]) == ("ENST11", "canonical", "0")
    assert r["acmg_label"] == "Likely_pathogenic" and r["is_acmg_PLP"] == "1"
    assert r["acmg_label_mane"] == "" and r["mane_discordant"] == ""
    # VEP-path order also has canonical before tsl.
    v, _, _ = run(tmp_path, [rec], "--pick-order", P.VEP_PATH_PICK_ORDER)
    assert v[0]["transcript_id"] == "ENST11"


def test_two_overlapping_genes_each_get_a_row(tmp_path):
    rec = ("17", 300, [
        blk("missense_variant", "GA", "GA_ID", "ENSTA", "LP", "PM1", mane="NM_A"),
        blk("intron_variant", "GB", "GB_ID", "ENSTB1", "VUS", "", canonical="YES"),
        blk("synonymous_variant", "GB", "GB_ID", "ENSTB2", "LB", "BP7"),
    ], ("0/1", "1/1", "0/0"))
    v, c, err = run(tmp_path, [rec])
    by = {r["gene"]: r for r in v}
    assert set(by) == {"GA", "GB"}
    assert by["GA"]["pick_tier"] == "only_candidate" and by["GA"]["is_mane_select"] == "1"
    assert by["GB"]["transcript_id"] == "ENSTB1" and by["GB"]["pick_tier"] == "canonical"
    assert {(x["gene"], x["person_id"], x["zygosity"]) for x in c} == {("GA", "S1", "het"),
                                                                       ("GA", "S2", "hom")}
    assert "1 variants with 2+ gene rows" in err


def test_flanking_only_neighbour_gene_dropped(tmp_path):
    rec = ("2", 400, [
        blk("missense_variant", "GMAIN", "GM", "ENSTM", "LP", "PM1", mane="NM_M"),
        blk("upstream_gene_variant", "GNEIGH", "GN", "ENSTN1", "VUS", "", mane="NM_N",
            impact="MODIFIER"),
        blk("downstream_gene_variant", "GNEIGH", "GN", "ENSTN2", "VUS", "", impact="MODIFIER"),
    ], HOMREF)
    v, _, err = run(tmp_path, [rec])
    assert [r["gene"] for r in v] == ["GMAIN"]
    assert "flanking-only gene rows dropped: 1 " in err
    v, _, err = run(tmp_path, [rec], "--keep-flanking-genes")
    assert sorted(r["gene"] for r in v) == ["GMAIN", "GNEIGH"]
    assert "flanking-only gene rows dropped: 0 " in err


def test_flanking_gene_with_any_other_consequence_is_kept(tmp_path):
    rec = ("2", 410, [
        blk("upstream_gene_variant", "G", "G", "ENST1", "VUS", "", impact="MODIFIER"),
        blk("intron_variant", "G", "G", "ENST2", "VUS", "", impact="MODIFIER"),
    ], HOMREF)
    v, _, _ = run(tmp_path, [rec])
    assert [r["transcript_id"] for r in v] == ["ENST2"]   # rank: intron(28) < upstream(32)


def test_same_gene_tie_broken_by_transcript_id(tmp_path):
    rec = ("3", 500, [
        blk("missense_variant", "GT", "GT", "ENST0009", "VUS", ""),
        blk("missense_variant", "GT", "GT", "ENST0002", "LP", "PM1"),
        blk("missense_variant", "GT", "GT", "ENST0005", "VUS", ""),
    ], HOMREF)
    v, _, err = run(tmp_path, [rec])
    assert v[0]["transcript_id"] == "ENST0002" and v[0]["pick_tier"] == "tie_id"
    # Order of blocks in the CSQ must not matter.
    rec2 = (rec[0], rec[1], list(reversed(rec[2])), rec[3])
    v2, _, _ = run(tmp_path, [rec2])
    assert v2[0]["transcript_id"] == "ENST0002"
    assert "tie_id=1" in err


def test_no_gene_variant_gets_sentinel(tmp_path):
    rec = ("4", 600, ["T|intergenic_variant|MODIFIER||||||||||||||||"], HOMREF)
    v, _, err = run(tmp_path, [rec])
    assert len(v) == 1 and v[0]["acmg_label"] == P.NO_TRANSCRIPT_SENTINEL
    assert v[0]["consequence"] == "intergenic_variant"
    assert "NO_TRANSCRIPT (no gene at all): 1" in err


def test_haploid_gt_is_hemi(tmp_path):
    rec = ("X", 700, [blk("stop_gained", "DMD", "GD", "ENSTD", "P", "PVS1", mane="NM_D",
                          impact="HIGH")], ("1", "0", "0/1"))
    _, c, err = run(tmp_path, [rec])
    assert {(x["person_id"], x["GT"], x["zygosity"]) for x in c} == {("S1", "1", "hemi"),
                                                                     ("S3", "0/1", "het")}
    assert "hemizygous calls: 1" in err


@pytest.mark.parametrize("gt,expected", [
    ("1", "hemi"), ("2", "hemi"), ("0", None), (".", None),
    ("0/1", "het"), ("1|0", "het"), ("1/1", "hom"), ("0/0", None), ("./.", None), ("./1", None),
])
def test_zygosity(gt, expected):
    assert _zygosity(gt) == expected


def test_tier_scoring_matches_pick_rs():
    assert P.parse_tsl("1") == 1 and P.parse_tsl("") == 255 and P.parse_tsl("tsl3") == 3
    assert P.parse_tsl("2 (assigned to previous version 1)") == 2
    assert P.consequence_rank("splice_region_variant&intron_variant") == 15
    assert P.consequence_rank("") == 2**32 - 1
    assert P.consequence_rank("made_up_term") < 2**32 - 1
    assert P.appris_rank("P1") < P.appris_rank("principal") < P.appris_rank("ALT1") \
        < P.appris_rank("alternative") < P.appris_rank("weird") < P.appris_rank("")
    assert P.appris_rank("appris_principal_1") == P.appris_rank("P1") == 1


def test_parse_pick_order():
    assert P.parse_pick_order(P.DEFAULT_PICK_ORDER) == [
        "mane_select", "mane_plus_clinical", "canonical", "tsl", "biotype", "ccds", "rank"]
    assert P.parse_pick_order("MANE, rank") == ["mane_select", "rank"]
    for bad in ("length", "rank,rank", "foo", ""):
        with pytest.raises(ValueError):
            P.parse_pick_order(bad)
