import csv
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin" / "build_carrier_matrix.py"


def _write(p, text):
    p.write_text(text, encoding="utf-8")


def test_gt_and_zygosity_roundtrip(tmp_path):
    # One qualifying ClinVar P/LP variant on chr10, and one on chrX for hemizygous.
    clinvar = tmp_path / "c.clinvar_plp.tsv"
    _write(clinvar,
        "chr\tpos\tref\talt\tgene\tclnsig\tclnrevstat\tcondition\tclnsigconf\tstars\tis_clinvar_PLP\n"
        "10\t100\tA\tG\tGENE1\tPathogenic\treviewed_by_expert_panel\tcond\t\t3\t1\n"
        "X\t200\tC\tT\tGENE2\tPathogenic\treviewed_by_expert_panel\tcond\t\t3\t1\n")

    gt = tmp_path / "gt.tsv"
    _write(gt,
        "10\t100\tA\tG\tS_HET\t0/1\n"
        "10\t100\tA\tG\tS_HOM\t1/1\n"
        "10\t100\tA\tG\tS_REF\t0/0\n"        # dropped (not a carrier)
        "X\t200\tC\tT\tS_HEMI\t1\n")          # haploid → hemizygous

    out = tmp_path / "carrier_matrix.tsv"
    subprocess.run(
        [sys.executable, str(BIN), "--gt", str(gt), "--clinvar", str(clinvar),
         "--out-long", str(out)],
        check=True,
    )

    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    by = {(r["person_id"]): r for r in rows}
    assert "S_REF" not in by                       # hom-ref dropped
    assert by["S_HET"]["GT"] == "0/1" and by["S_HET"]["zygosity"] == "het"
    assert by["S_HOM"]["GT"] == "1/1" and by["S_HOM"]["zygosity"] == "hom_alt"
    assert by["S_HEMI"]["GT"] == "1" and by["S_HEMI"]["zygosity"] == "hemizygous"
    # header contract
    assert list(rows[0].keys())[:8] == [
        "chr", "pos", "ref", "alt", "gene", "person_id", "GT", "zygosity"]


def test_zygosity_helper():
    sys.path.insert(0, str(BIN.parent))
    import build_carrier_matrix as bcm
    assert bcm.zygosity("0/1", "10") == "het"
    assert bcm.zygosity("1|1", "10") == "hom_alt"
    assert bcm.zygosity("1", "X") == "hemizygous"
    assert bcm.zygosity("1/.", "X") == "hemizygous"
    assert bcm.zygosity("./.", "10") == "missing"
