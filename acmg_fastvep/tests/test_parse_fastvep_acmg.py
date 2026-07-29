import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARSER = HERE.parent / "parse_fastvep_acmg.py"

# fastVEP CSQ header: SYMBOL ... ACMG, ACMG_CRITERIA are the relevant fields.
CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from fastVEP. '
    'Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|BIOTYPE|ACMG|ACMG_CRITERIA">'
)


def _vcf(tmp_path):
    p = tmp_path / "fv.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        CSQ_HEADER,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        # BRCA1 frameshift: two transcript blocks, LP is the most severe.
        "17\t43000000\t.\tAC\tA\t.\t.\tCSQ=A|frameshift_variant|HIGH|BRCA1|ENSG|ENST1|protein_coding|LP|PVS1&PM2_Supporting,"
        "A|frameshift_variant|HIGH|BRCA1|ENSG|ENST2|protein_coding|VUS|PM2_Supporting",
        # TP53 missense: Pathogenic.
        "17\t7670000\t.\tG\tA\t.\t.\tCSQ=A|missense_variant|MODERATE|TP53|ENSG|ENST3|protein_coding|P|PS1&PM1&PP3_Moderate",
        # Benign standalone.
        "1\t100\t.\tC\tT\t.\t.\tCSQ=T|missense_variant|MODERATE|GENEB|ENSG|ENST4|protein_coding|B|BA1",
        # No ACMG call on any block → skipped.
        "2\t200\t.\tG\tC\t.\t.\tCSQ=C|intron_variant|MODIFIER|GENEI|ENSG|ENST5|protein_coding||",
    ]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def test_parse_fastvep_acmg(tmp_path):
    vcf = _vcf(tmp_path)
    out = tmp_path / "acmg_plp.tsv"
    subprocess.run([sys.executable, str(PARSER), "--vcf", str(vcf), "--out", str(out)], check=True)
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    by = {(r["chr"], r["pos"]): r for r in rows}

    # header/schema matches the pipeline's acmg_plp.tsv
    assert list(rows[0].keys()) == [
        "chr", "pos", "ref", "alt", "gene", "acmg_label", "acmg_criteria",
        "n_pathogenic_criteria", "n_benign_criteria", "is_acmg_PLP"]

    # BRCA1: most-severe across transcripts = LP (not VUS)
    r = by[("17", "43000000")]
    assert r["gene"] == "BRCA1" and r["acmg_label"] == "Likely_pathogenic"
    assert r["acmg_criteria"] == "PVS1;PM2_Supporting"
    assert r["n_pathogenic_criteria"] == "2" and r["n_benign_criteria"] == "0"
    assert r["is_acmg_PLP"] == "1"

    # TP53: Pathogenic
    r = by[("17", "7670000")]
    assert r["acmg_label"] == "Pathogenic" and r["is_acmg_PLP"] == "1"
    assert r["n_pathogenic_criteria"] == "3"

    # Benign: BA1, not P/LP
    r = by[("1", "100")]
    assert r["acmg_label"] == "Benign" and r["is_acmg_PLP"] == "0"
    assert r["n_benign_criteria"] == "1"

    # variant with no ACMG call is dropped
    assert ("2", "200") not in by
