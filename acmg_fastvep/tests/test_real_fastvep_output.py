"""The carrier builders against REAL fastVEP v0.5.0 output.

Every other test in this directory feeds the builders CSQ lines typed by hand.
These run them on what the released fastVEP v0.5.0 binary wrote for a synthetic
cohort (tests/synthetic/real_fastvep_v050/README.md says how it was made), and compare with a hand count that
is written down in this file from how the cohort was designed, not from what the
builder printed.

What this does and does not exercise: the real CSQ header and block layout of
fastVEP v0.5.0 (51 fields, `PICK`, `ACMG`, `ACMG_CRITERIA`, `ACMG_DISORDERS`),
real consequence terms, real PICK flags and real ACMG labels from a three-source
stack (ClinGen-style gene validity, gene constraint, one gnomAD record). It does
not exercise a real GFF3, ClinVar, REVEL or SpliceAI, and the genes are synthetic.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PICK = HERE.parent / "build_acmg_carriers_pick.py"
MANE = HERE.parent / "build_acmg_carriers.py"
# under tests/synthetic/ because the repository ignores *.vcf and hooks/pre-commit exempts only that tree
DATA = HERE.parents[1] / "tests" / "synthetic" / "real_fastvep_v050"
COHORT = DATA / "real_fastvep_v0.5.0_cohort.vcf"
SPLIT = DATA / "real_fastvep_v0.5.0_multiallelic_split.vcf"

# Ten samples S1..S10. At the two "zoo" sites they carry, in this order:
#   0/1  1/1  1  ./1  0|1  ./.  0/0  1/.  1|1  0
# Definite carriers: S1 het, S2 hom, S3 hemi (haploid alt), S5 het (phased), S9 hom (phased).
# Not carriers by this pipeline's rule: S4 ./1 and S8 1/. (half calls), S6 ./., S7 0/0, S10 haploid ref.
ZOO = {"S1": "het", "S2": "hom", "S3": "hemi", "S5": "het", "S9": "hom"}
ONE = {"S1": "het"}
# (chr, pos, gene) -> {person: zygosity}, from the design of the cohort:
TRUTH = {
    ("1", "1054", "GENEA"): ZOO,   # stop_gained in a gene with MANE Select
    ("1", "1100", "GENEA"): ONE,   # frameshift insertion, same gene
    ("1", "5054", "GENEC"): ONE,   # stop_gained, gene WITHOUT a MANE transcript, two transcripts: canonical is picked
    ("1", "8054", "GENED"): ONE,   # stop_gained, MANE Plus Clinical transcript only
    ("1", "20054", "GENEH"): ONE,  # stop_gained, isolated gene without MANE
    ("1", "25054", "GENEI"): ONE,  # stop_gained, gene without MANE beside a MANE gene
    ("X", "1054", "GENEX"): ZOO,   # same zoo on chrX
}
# No carriers by design: 1:1160 (missense in GENEA, 5' UTR of GENEB, both VUS) and 1:12304 (stop_gained only
# in a non-MANE transcript of GENEF whose MANE transcript has an intron variant there: VUS under a pick).
# (pos, gene) -> (is_acmg_PLP, acmg_label_mane) for the gene rows of the cohort without its multi-allelic site:
GENE_ROWS = {
    ("1054", "GENEA"): ("1", "Likely_pathogenic"), ("1100", "GENEA"): ("1", "Likely_pathogenic"),
    ("1160", "GENEA"): ("0", "Uncertain_significance"), ("1160", "GENEB"): ("0", "Uncertain_significance"),
    ("5054", "GENEC"): ("1", "NO_MANE_TRANSCRIPT"), ("8054", "GENED"): ("1", "NO_MANE_TRANSCRIPT"),
    ("12304", "GENEF"): ("0", "Uncertain_significance"),
    ("20054", "GENEH"): ("1", "NO_MANE_TRANSCRIPT"), ("25054", "GENEI"): ("1", "NO_MANE_TRANSCRIPT"),
    ("1054", "GENEX"): ("1", "Likely_pathogenic"),
}


def _run(script, vcf, tmp_path, *extra):
    ov, oc = tmp_path / "variants.tsv", tmp_path / "carriers.tsv"
    cmd = [sys.executable, str(script), "--vcf", str(vcf), "--out-variants", str(ov), "--out-carriers", str(oc), *extra]
    if script == PICK:
        cmd.append("--no-markers")
    res = subprocess.run(cmd, capture_output=True, text=True)
    return res, ov, oc


def _biallelic_cohort(tmp_path):
    """The cohort fixture without its one multi-allelic record (1:16054 C>T,G)."""
    out = tmp_path / "cohort_biallelic.vcf"
    keep = [ln for ln in COHORT.read_text(encoding="utf-8").splitlines(keepends=True)
            if ln.startswith("#") or ln.split("\t")[4].count(",") == 0]
    out.write_text("".join(keep), encoding="utf-8")
    return out


def _rows(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _carriers(path):
    out = {}
    for r in _rows(path):
        out.setdefault((r["chr"], r["pos"], r["gene"]), {})[r["person_id"]] = r["zygosity"]
    return out


def test_pick_carriers_equal_the_hand_count_on_real_output(tmp_path):
    res, ov, oc = _run(PICK, _biallelic_cohort(tmp_path), tmp_path)
    assert res.returncode == 0, res.stderr
    assert _carriers(oc) == TRUTH


def test_pick_gene_rows_on_real_output(tmp_path):
    res, ov, oc = _run(PICK, _biallelic_cohort(tmp_path), tmp_path)
    assert res.returncode == 0, res.stderr
    got = {(r["pos"], r["gene"]): (r["is_acmg_PLP"], r["acmg_label_mane"]) for r in _rows(ov) if r["chr"] in ("1", "X")}
    assert got == GENE_ROWS
    # every picked row is fastVEP's own flagged block, and nothing was flagged as an anomaly
    assert {r["pick_tier"] for r in _rows(ov)} == {"fastvep_pick"}
    assert {r["pick_flag"] for r in _rows(ov)} == {"ok"}


def test_pick_builder_refuses_a_multiallelic_record_by_default(tmp_path):
    res, ov, oc = _run(PICK, COHORT, tmp_path)
    assert res.returncode == 2, res.stderr
    assert "1:16054" in res.stderr and "bcftools norm -m -any" in res.stderr


def test_pick_builder_skips_a_multiallelic_record_when_told_to_and_guesses_nothing(tmp_path):
    res, ov, oc = _run(PICK, COHORT, tmp_path, "--allow-multiallelic")
    assert res.returncode == 0, res.stderr
    assert "multi-allelic ALT records skipped: 1" in res.stderr
    assert all(r["pos"] != "16054" for r in _rows(ov)) and all(r["pos"] != "16054" for r in _rows(oc))


def test_split_form_of_the_multiallelic_site_gives_exactly_the_true_carriers(tmp_path):
    # `bcftools norm -m -any` (the pipeline's NORM_QC) turned C>T,G with samples 0/2 0/1 1/2 2/2 into two records:
    #   C>T: 0/0 0/1 1/0 0/0   and   C>G: 0/1 0/0 0/1 1/1.
    # Carriers of the stop_gained allele (T) are S2 and S3, and only them; the G allele is a VUS and has none.
    res, ov, oc = _run(PICK, SPLIT, tmp_path)
    assert res.returncode == 0, res.stderr
    assert _carriers(oc) == {("1", "16054", "GENEG"): {"S2": "het", "S3": "het"}}
    plp = {(r["alt"], r["is_acmg_PLP"]) for r in _rows(ov) if r["gene"] == "GENEG"}
    assert plp == {("T", "1"), ("G", "0")}


def test_mane_only_builder_refuses_a_multiallelic_record_by_default(tmp_path):
    res, ov, oc = _run(MANE, COHORT, tmp_path)
    assert res.returncode == 2, res.stderr
    assert "1:16054" in res.stderr and "bcftools norm -m -any" in res.stderr


def test_mane_only_builder_skips_it_when_told_to(tmp_path):
    res, ov, oc = _run(MANE, COHORT, tmp_path, "--allow-multiallelic")
    assert res.returncode == 0, res.stderr
    assert "multi-allelic ALT records skipped: 1" in res.stderr
    assert all(r["pos"] != "16054" for r in _rows(ov)) and all(r["pos"] != "16054" for r in _rows(oc))
    # the zoo sites are unaffected
    got = _carriers(oc)
    assert got[("1", "1054", "GENEA")] == ZOO and got[("X", "1054", "GENEX")] == ZOO
