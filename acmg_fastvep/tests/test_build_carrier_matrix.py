"""Functional test for build_carrier_matrix.sh (standalone carrier assembly).

Exercises the whole orchestration without bcftools (falls back to
fixture_gt_extract.py): a variant that's ClinVar-P/LP and a *different* variant
that's ACMG-P/LP must BOTH yield carrier rows, with the right flags — proving
ACMG-only variants (absent from a ClinVar-only matrix) are picked up.
"""
from __future__ import annotations
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "build_carrier_matrix.sh"

VCF = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\t1/1
1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/0\t0/1\t0/0
1\t300\t.\tG\tA\t.\t.\t.\tGT\t0/0\t0/0\t0/0
"""

# ClinVar P/LP at pos 100; ACMG P/LP at a DIFFERENT variant, pos 200.
CLINVAR = (
    "chr\tpos\tref\talt\tgene\tclnsig\tclnrevstat\tcondition\tclnsigconf\tstars\tis_clinvar_PLP\n"
    "1\t100\tA\tG\tGENEA\tPathogenic\treviewed\t\t\t3\t1\n"
    "1\t200\tC\tT\tGENEB\t\t\t\t\t0\t0\n"
)
ACMG = (
    "chr\tpos\tref\talt\tgene\tacmg_label\tacmg_criteria\tn_pathogenic_criteria\tn_benign_criteria\tis_acmg_PLP\n"
    "1\t200\tC\tT\tGENEB\tLikely_pathogenic\tPVS1;PM2_Supporting\t2\t0\t1\n"
    "1\t100\tA\tG\tGENEA\tUncertain_significance\t\t0\t0\t0\n"
)


def test_combined_carrier_matrix(tmp_path):
    norm = tmp_path / "norm_qc"; norm.mkdir()
    cvd = tmp_path / "variants"; cvd.mkdir()
    acd = tmp_path / "fastvep"; acd.mkdir()
    (norm / "chunk.norm.vcf").write_text(VCF)          # uncompressed → python fallback
    (cvd / "chunk.clinvar_plp.tsv").write_text(CLINVAR)
    (acd / "chunk.acmg_plp.tsv").write_text(ACMG)
    out = tmp_path / "carrier_matrix.tsv"

    r = subprocess.run(
        ["bash", str(SCRIPT), "--norm-dir", str(norm), "--clinvar-dir", str(cvd),
         "--acmg-dir", str(acd), "-o", str(out), "--bcftools", "definitely_not_bcftools"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, f"script failed:\n{r.stdout}\n{r.stderr}"
    rows = [ln.split("\t") for ln in out.read_text().splitlines()]
    hdr = rows[0]
    assert hdr[:6] == ["chr", "pos", "ref", "alt", "gene", "person_id"]
    ci = {name: i for i, name in enumerate(hdr)}
    data = rows[1:]

    # pos 100 (ClinVar P/LP): carriers S1 (0/1) and S3 (1/1), flagged clinvar not acmg
    p100 = [r for r in data if r[1] == "100"]
    assert {r[ci["person_id"]] for r in p100} == {"S1", "S3"}
    assert all(r[ci["is_clinvar_PLP"]] == "1" and r[ci["is_acmg_PLP"]] == "0" for r in p100)

    # pos 200 (ACMG-ONLY P/LP): carrier S2 (0/1), flagged acmg not clinvar — the key case
    p200 = [r for r in data if r[1] == "200"]
    assert {r[ci["person_id"]] for r in p200} == {"S2"}
    assert all(r[ci["is_acmg_PLP"]] == "1" and r[ci["is_clinvar_PLP"]] == "0" for r in p200)

    # pos 300 (no P/LP, all hom-ref): no rows
    assert not [r for r in data if r[1] == "300"]
