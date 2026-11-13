"""build_acmg_carriers_pick.py --selection fastvep_pick (the default) on fastVEP
v0.5.0-shaped output: CSQ with PICK after FLAGS (--flag-pick-allele-gene) and
ACMG_DISORDERS as the last field (--acmg), as crates/fastvep-io/src/output.rs
lays them out at the v0.5.0 tag."""
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "build_acmg_carriers_pick.py"
DRIVER = HERE.parent / "run_acmg_carriers_chr.sh"

FIELDS = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
          "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "FLAGS", "PICK", "CANONICAL",
          "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "CCDS", "ACMG", "ACMG_CRITERIA",
          "ACMG_DISORDERS"]


def header(fields=FIELDS):
    return ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from '
            f'fastVEP. Format: {"|".join(fields)}">')


def blk(cons, sym, gene, feat, acmg, crit="", *, pick="", mane="", canonical="", tsl="",
        disorders="", impact="MODERATE", allele="T", fields=FIELDS):
    v = {"Allele": allele, "Consequence": cons, "IMPACT": impact, "SYMBOL": sym, "Gene": gene,
         "Feature_type": "Transcript", "Feature": feat, "BIOTYPE": "protein_coding",
         "PICK": pick, "CANONICAL": canonical, "MANE_SELECT": mane, "TSL": tsl,
         "ACMG": acmg, "ACMG_CRITERIA": crit, "ACMG_DISORDERS": disorders}
    return "|".join(v.get(f, "") for f in fields)


def write_vcf(path: Path, records, fields=FIELDS, gz=False):
    lines = ["##fileformat=VCFv4.2", header(fields),
             "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                        "S1", "S2"])]
    for chrom, pos, blocks, gts in records:
        lines.append("\t".join([chrom, str(pos), ".", "C", "T", ".", "PASS",
                                "CSQ=" + ",".join(blocks), "GT", *gts]))
    text = "\n".join(lines) + "\n"
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def run(tmp_path, records, *extra, fields=FIELDS, check=True, rc=0):
    vcf = write_vcf(tmp_path / "in.fastvep.vcf", records, fields)
    ov, oc = tmp_path / "v.tsv", tmp_path / "c.tsv"
    res = subprocess.run([sys.executable, str(SCRIPT), "--vcf", str(vcf), "--out-variants",
                          str(ov), "--out-carriers", str(oc), *extra],
                         capture_output=True, text=True)
    if check:
        assert res.returncode == rc, res.stderr
        return (list(csv.DictReader(ov.open(), delimiter="\t")),
                list(csv.DictReader(oc.open(), delimiter="\t")), res.stderr)
    return res


def test_takes_fastveps_pick_not_a_python_pick(tmp_path):
    """fastVEP flagged the non-MANE stop_gained (e.g. a rank-first --pick-order);
    the row must report that block, not the MANE one a Python port would choose."""
    rec = ("17", 100, [
        blk("missense_variant", "BRCA1", "G1", "ENST_MANE", "VUS", "PM2", mane="NM_1.1",
            canonical="YES", tsl="1", disorders="AD:VUS:PM2"),
        blk("stop_gained", "BRCA1", "G1", "ENST_OTHER", "P", "PVS1&PM2", pick="1", tsl="1",
            impact="HIGH", disorders="AD:P:PVS1+PM2&AR:LP:PVS1"),
    ], ("0/1", "0/0"))
    v, c, err = run(tmp_path, [rec])
    assert len(v) == 1
    r = v[0]
    assert (r["transcript_id"], r["pick_tier"], r["pick_flag"]) == ("ENST_OTHER", "fastvep_pick", "ok")
    assert r["acmg_label"] == "Pathogenic" and r["is_acmg_PLP"] == "1"
    assert r["acmg_criteria"] == "PVS1;PM2"
    assert r["acmg_disorders"] == "AD:P:PVS1+PM2&AR:LP:PVS1"
    assert r["is_mane_select"] == "0"
    assert r["acmg_label_mane"] == "Uncertain_significance" and r["mane_discordant"] == "1"
    assert [x["person_id"] for x in c] == ["S1"]
    assert "selection=fastvep_pick" in err and "DEPRECATED" not in err


def test_no_mane_gene_says_so(tmp_path):
    rec = ("1", 200, [blk("stop_gained", "GX", "GX_ID", "ENST1", "LP", "PVS1", pick="1",
                          canonical="YES", impact="HIGH")], ("0/0", "1/1"))
    v, c, _ = run(tmp_path, [rec])
    assert v[0]["acmg_label_mane"] == "NO_MANE_TRANSCRIPT" and v[0]["mane_discordant"] == ""
    assert [(x["person_id"], x["zygosity"]) for x in c] == [("S2", "hom")]


def test_zero_or_two_picks_are_flagged_not_resolved(tmp_path):
    none = ("2", 300, [blk("stop_gained", "GA", "GA_ID", "ENSTA1", "P", "PVS1", mane="NM_A"),
                       blk("missense_variant", "GA", "GA_ID", "ENSTA2", "VUS")], ("0/1", "0/1"))
    two = ("2", 400, [blk("stop_gained", "GB", "GB_ID", "ENSTB2", "P", "PVS1", pick="1"),
                      blk("missense_variant", "GB", "GB_ID", "ENSTB1", "VUS", pick="1")],
           ("0/1", "0/0"))
    # Outputs are written, and the script then fails: per pick.rs this cannot happen.
    v, c, err = run(tmp_path, [none, two], rc=3)
    assert "ERROR: 2 (allele, gene) row(s) had no or several PICK=1 blocks" in err
    by = {r["gene"]: r for r in v}
    assert by["GA"]["acmg_label"] == "NO_PICK_FLAG" and by["GA"]["is_acmg_PLP"] == "0"
    assert by["GA"]["pick_flag"] == "no_pick" and by["GA"]["transcript_id"] == ""
    assert by["GA"]["acmg_label_mane"] == "Pathogenic"      # still shown, for review
    assert by["GB"]["acmg_label"] == "MULTIPLE_PICK_FLAGS" and by["GB"]["is_acmg_PLP"] == "0"
    assert by["GB"]["pick_flag"] == "multi_pick:ENSTB1=Uncertain_significance,ENSTB2=Pathogenic"
    assert by["GB"]["transcript_id"] == "ENSTB2;ENSTB1"
    assert c == []                                          # flagged rows carry no carriers
    assert "with no PICK=1: 1; with more than one: 1" in err
    assert err.count("WARN:") == 2
    v2, _, err2 = run(tmp_path, [none, two], "--allow-pick-anomalies")
    assert v2 == v and "accepted (--allow-pick-anomalies)" in err2


def test_one_row_per_gene(tmp_path):
    rec = ("17", 500, [
        blk("missense_variant", "GA", "GA_ID", "ENSTA", "LP", "PM1", pick="1", mane="NM_A"),
        blk("intron_variant", "GB", "GB_ID", "ENSTB1", "VUS", pick="1", canonical="YES"),
        blk("synonymous_variant", "GB", "GB_ID", "ENSTB2", "LB", "BP7"),
    ], ("0/1", "0/0"))
    v, c, _ = run(tmp_path, [rec])
    assert sorted((r["gene"], r["transcript_id"]) for r in v) == [("GA", "ENSTA"), ("GB", "ENSTB1")]
    assert [(x["gene"], x["person_id"]) for x in c] == [("GA", "S1")]


def test_missing_pick_field_is_an_error(tmp_path):
    fields = [f for f in FIELDS if f != "PICK"]
    rec = ("1", 1, [blk("stop_gained", "G", "G", "E", "P", fields=fields)], ("0/1", "0/0"))
    res = run(tmp_path, [rec], fields=fields, check=False)
    assert res.returncode == 2 and "no PICK field" in res.stderr
    # the deprecated port still reads it, and says it is deprecated
    v, _, err = run(tmp_path, [rec], "--selection", "python_pick", fields=fields)
    assert v[0]["acmg_label"] == "Pathogenic" and "DEPRECATED" in err


def test_pick_order_rejected_with_fastvep_pick(tmp_path):
    rec = ("1", 1, [blk("stop_gained", "G", "G", "E", "P", pick="1")], ("0/1", "0/0"))
    res = run(tmp_path, [rec], "--pick-order", "rank", check=False)
    assert res.returncode == 2 and "python_pick only" in res.stderr


def test_chr_driver_defaults_to_fastvep_pick(tmp_path):
    indir = tmp_path / "results-fastvep-chr21_fv050"
    indir.mkdir()
    write_vcf(indir / "c1.fastvep.vcf.gz", [("21", 10, [
        blk("stop_gained", "G1", "G1", "E1", "P", "PVS1", pick="1", impact="HIGH")], ("0/1", "0/0"))],
        gz=True)
    write_vcf(indir / "c2.fastvep.vcf.gz", [("21", 20, [
        blk("missense_variant", "G2", "G2", "E2", "VUS"),
        blk("missense_variant", "G2", "G2", "E3", "LP", pick="1", mane="NM_3")], ("0/0", "0/1"))],
        gz=True)
    out = tmp_path / "results-carriers-fastvep_pick-chr21_fv050"
    res = subprocess.run(["bash", str(DRIVER), "--in-dir", str(indir), "-o", str(out), "-j", "2"],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    v = list(csv.DictReader((out / "acmg_variants.tsv").open(), delimiter="\t"))
    assert sorted((r["gene"], r["transcript_id"], r["pick_tier"]) for r in v) == \
        [("G1", "E1", "fastvep_pick"), ("G2", "E3", "fastvep_pick")]
    summary = (out / "summary.txt").read_text()
    assert "selection=fastvep_pick" in summary and "2 variants -> 2 (variant, gene) rows" in summary
    bad = subprocess.run(["bash", str(DRIVER), "--in-dir", str(indir), "-o", str(tmp_path / "x"),
                          "--pick-order", "rank"], capture_output=True, text=True)
    assert bad.returncode == 2


def test_compare_reads_fastvep_pick_dirs_across_suffixes(tmp_path):
    """MANE (v6 run) vs fastvep_pick (fv050 run) on chr21: dirs differ in prefix
    and suffix, and are matched with --pick-prefix / --strip-dir-suffix."""
    vcf = write_vcf(tmp_path / "in.fastvep.vcf", [("21", 10, [
        blk("missense_variant", "G1", "G1", "EM", "VUS", mane="NM_1"),
        blk("stop_gained", "G1", "G1", "EO", "P", "PVS1", pick="1", impact="HIGH")],
        ("0/1", "0/0"))])
    md = tmp_path / "mane" / "results-carriers-chr21_v6"
    pd = tmp_path / "fv" / "results-carriers-fastvep_pick-chr21_fv050"
    for script, d in (("build_acmg_carriers.py", md), ("build_acmg_carriers_pick.py", pd)):
        d.mkdir(parents=True)
        subprocess.run([sys.executable, str(HERE.parent / script), "--vcf", str(vcf),
                        "--out-variants", str(d / "acmg_variants.tsv"),
                        "--out-carriers", str(d / "acmg_carriers.tsv")], check=True,
                       capture_output=True)
    out = tmp_path / "fv" / "comparison"
    res = subprocess.run([sys.executable, str(HERE.parent / "compare_acmg_pick_vs_mane.py"),
                          "--mane-root", str(tmp_path / "mane"), "--pick-root", str(tmp_path / "fv"),
                          "--pick-prefix", "results-carriers-fastvep_pick-", "--strip-dir-suffix",
                          "--out", str(out)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert "comparing 1 chromosome(s): ['chr21']" in res.stderr
    tiers = (out / "pick_tier_counts.tsv").read_text()
    assert "fastvep_pick" in tiers
    assert "EO" in (out / "mane_discordant_plp.tsv").read_text()


def test_chr_driver_fails_on_pick_anomalies_after_writing(tmp_path):
    indir = tmp_path / "results-fastvep-chr21_fv050"
    indir.mkdir()
    write_vcf(indir / "c1.fastvep.vcf.gz", [("21", 10, [
        blk("stop_gained", "G1", "G1", "E1", "P", "PVS1", pick="1", impact="HIGH")], ("0/1", "0/0"))],
        gz=True)
    write_vcf(indir / "c2.fastvep.vcf.gz", [("21", 20, [
        blk("missense_variant", "G2", "G2", "E2", "VUS"),
        blk("missense_variant", "G2", "G2", "E3", "LP", mane="NM_3")], ("0/0", "0/1"))], gz=True)
    out = tmp_path / "out"
    res = subprocess.run(["bash", str(DRIVER), "--in-dir", str(indir), "-o", str(out), "-j", "2"],
                         capture_output=True, text=True)
    assert res.returncode == 5, res.stderr
    assert "without exactly one PICK=1" in res.stderr
    failed = (out / "FAILED").read_text()
    assert "exit_code\t5" in failed and "PICK anomalies" in failed and "no PICK=1: 1" in failed
    assert not (out / "DONE").exists()
    v = list(csv.DictReader((out / "acmg_variants.tsv").open(), delimiter="\t"))
    assert sorted(r["acmg_label"] for r in v) == ["NO_PICK_FLAG", "Pathogenic"]
    assert "with no PICK=1: 1" in (out / "summary.txt").read_text()
    ok = subprocess.run(["bash", str(DRIVER), "--in-dir", str(indir), "-o", str(out), "-j", "2",
                         "--allow-pick-anomalies"], capture_output=True, text=True)
    assert ok.returncode == 0, ok.stderr
    assert (out / "DONE").exists() and not (out / "FAILED").exists()
    assert "allow_pick_anomalies=1" in (out / "DONE").read_text()


def _carrier_dir(root: Path, name: str, records, *extra):
    d = root / name
    d.mkdir(parents=True)
    vcf = write_vcf(d.parent / f"{name}.vcf", records)
    return d, subprocess.run([sys.executable, str(SCRIPT), "--vcf", str(vcf), "--out-variants",
                              str(d / "acmg_variants.tsv"), "--out-carriers",
                              str(d / "acmg_carriers.tsv"), *extra], capture_output=True, text=True)


GOOD = [("21", 10, [blk("stop_gained", "BRCA2", "G1", "E1", "P", "PVS1", pick="1", mane="NM_1",
                        impact="HIGH")], ("0/1", "0/0"))]
BAD = [("21", 20, [blk("stop_gained", "BRCA2", "G1", "E1", "P", "PVS1", mane="NM_1"),
                   blk("missense_variant", "BRCA2", "G1", "E2", "VUS")], ("0/1", "0/0"))]


def test_builder_writes_markers(tmp_path):
    d, r = _carrier_dir(tmp_path, "ok", GOOD)
    assert r.returncode == 0 and (d / "DONE").exists() and not (d / "FAILED").exists()
    d2, r2 = _carrier_dir(tmp_path, "bad", BAD)
    assert r2.returncode == 3 and (d2 / "FAILED").exists() and not (d2 / "DONE").exists()
    assert "exit_code\t3" in (d2 / "FAILED").read_text()
    # a later successful run in the same dir clears the stale FAILED
    vcf = write_vcf(tmp_path / "good.vcf", GOOD)
    r3 = subprocess.run([sys.executable, str(SCRIPT), "--vcf", str(vcf), "--out-variants",
                         str(d2 / "acmg_variants.tsv"), "--out-carriers", str(d2 / "acmg_carriers.tsv")],
                        capture_output=True, text=True)
    assert r3.returncode == 0 and (d2 / "DONE").exists() and not (d2 / "FAILED").exists()


def test_readers_refuse_failed_dirs(tmp_path):
    root = tmp_path / "fv"
    _carrier_dir(root, "results-carriers-fastvep_pick-chr21_fv050", BAD)
    genes = tmp_path / "genes.txt"
    genes.write_text("BRCA2\n")
    ex = [sys.executable, str(HERE.parent / "extract_gene_carriers.py"), "--pick-root", str(root),
          "--prefix", "results-carriers-fastvep_pick-", "--genes", str(genes),
          "--out", str(tmp_path / "o.tsv"), "--out-summary", str(tmp_path / "s.tsv")]
    r = subprocess.run(ex, capture_output=True, text=True)
    assert r.returncode == 4 and "FAILED marker" in r.stderr
    r = subprocess.run([*ex, "--allow-pick-anomalies"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    mane = tmp_path / "mane" / "results-carriers-chr21_v6"
    mane.mkdir(parents=True)
    vcf = write_vcf(tmp_path / "m.vcf", BAD)
    subprocess.run([sys.executable, str(HERE.parent / "build_acmg_carriers.py"), "--vcf", str(vcf),
                    "--out-variants", str(mane / "acmg_variants.tsv"),
                    "--out-carriers", str(mane / "acmg_carriers.tsv")], check=True, capture_output=True)
    cmp_ = [sys.executable, str(HERE.parent / "compare_acmg_pick_vs_mane.py"),
            "--mane-root", str(tmp_path / "mane"), "--pick-root", str(root),
            "--pick-prefix", "results-carriers-fastvep_pick-", "--strip-dir-suffix",
            "--out", str(root / "comparison")]
    r = subprocess.run(cmp_, capture_output=True, text=True)
    assert r.returncode == 4 and "FAILED marker" in r.stderr
    assert "no DONE marker" not in r.stderr or "mane" in r.stderr
    r = subprocess.run([*cmp_, "--allow-pick-anomalies"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_driver_writes_failed_on_missing_chunk(tmp_path, monkeypatch):
    """The existing exit-4 path (a chunk's output missing) also leaves FAILED."""
    indir = tmp_path / "results-fastvep-chr21_fv050"
    indir.mkdir()
    write_vcf(indir / "c1.fastvep.vcf.gz", GOOD, gz=True)
    (indir / "c2.fastvep.vcf.gz").write_bytes(b"")          # empty chunk -> its worker fails
    out = tmp_path / "out"
    res = subprocess.run(["bash", str(DRIVER), "--in-dir", str(indir), "-o", str(out), "-j", "1"],
                         capture_output=True, text=True)
    assert res.returncode != 0
    assert (out / "FAILED").exists() and not (out / "DONE").exists()
