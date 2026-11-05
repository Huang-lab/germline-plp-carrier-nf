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


def test_chr_driver_pick_mode(tmp_path):
    """run_acmg_carriers_chr.sh --mode pick: two gz chunks -> concatenated tables
    plus an aggregated summary.txt; refuses an outdir inside the input dir."""
    import gzip
    import shutil
    indir = tmp_path / "results-fastvep-chr17_v6"
    indir.mkdir()
    recs = [
        [("17", 100, [blk("stop_gained", "G1", "G1", "ENST1", "P", "PVS1", canonical="YES",
                          impact="HIGH")], ("0/1", "0/0", "0/0"))],
        [("17", 200, [blk("missense_variant", "G2", "G2", "ENST2", "VUS", "", mane="NM_2"),
                      blk("upstream_gene_variant", "G3", "G3", "ENST3", "VUS", "")], HOMREF)],
    ]
    for i, r in enumerate(recs):
        sub = tmp_path / f"c{i}"
        sub.mkdir()
        src = write_vcf(sub, r)
        with open(src, "rb") as fi, gzip.open(indir / f"chunk{i}.fastvep.vcf.gz", "wb") as fo:
            shutil.copyfileobj(fi, fo)
    driver = HERE.parent / "run_acmg_carriers_chr.sh"
    out = tmp_path / "results-acmg-pick" / "results-carriers-pick-chr17_v6"
    res = subprocess.run(["bash", str(driver), "--in-dir", str(indir), "-o", str(out),
                          "--mode", "pick", "-j", "2"], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    v = list(csv.DictReader((out / "acmg_variants.tsv").open(), delimiter="\t"))
    assert sorted((r["gene"], r["pick_tier"]) for r in v) == [("G1", "only_candidate"),
                                                              ("G2", "only_candidate")]
    summary = (out / "summary.txt").read_text()
    assert "2 variants -> 2 (variant, gene) rows" in summary
    assert "flanking-only gene rows dropped: 1 " in summary
    assert len((out / "acmg_carriers.tsv").read_text().splitlines()) == 2

    bad = subprocess.run(["bash", str(driver), "--in-dir", str(indir), "-o", str(indir / "x"),
                          "--mode", "pick"], capture_output=True, text=True)
    assert bad.returncode == 3 and "refusing" in bad.stderr
    assert not (indir / "x").exists()
    bad = subprocess.run(["bash", str(driver), "--in-dir", str(indir), "-o", str(out),
                          "--pick-order", "rank"], capture_output=True, text=True)
    assert bad.returncode == 2


def test_submit_pick_mode_with_stub_bsub(tmp_path):
    """submit_acmg_carriers.sh --mode pick -> run_acmg_carriers_task.lsf, with a
    stub bsub that runs every array index locally. Checks the output lands in
    --out-root/results-carriers-pick-<chr>/, the batch root is left untouched,
    and a --pick-order with commas survives the trip through bsub."""
    import gzip
    import os
    import shutil
    batch = tmp_path / "batch"
    for chrom, pos in (("chr21", 100), ("chr22", 200)):
        d = batch / f"results-fastvep-{chrom}_v6"
        d.mkdir(parents=True)
        src = write_vcf(tmp_path, [(chrom, pos, [
            blk("missense_variant", "G", "G", "ENST2", "VUS", "", canonical="YES"),
            blk("stop_gained", "G", "G", "ENST1", "P", "PVS1", impact="HIGH")],
            ("0/1", "0/0", "1"))])
        with open(src, "rb") as fi, gzip.open(d / "c.fastvep.vcf.gz", "wb") as fo:
            shutil.copyfileobj(fi, fo)
    before = sorted(str(p) for p in batch.rglob("*"))

    stub = tmp_path / "bin"
    stub.mkdir()
    (stub / "bsub").write_text(
        '#!/usr/bin/env bash\n'
        'script="$(cat)"\n'
        'n="$(printf "%s\\n" "$@" | sed -n "s/.*\\[1-\\([0-9]*\\)\\]$/\\1/p")"\n'
        'for i in $(seq 1 "$n"); do LSB_JOBINDEX=$i bash -c "$script" || exit 1; done\n')
    (stub / "bsub").chmod(0o755)
    out_root = tmp_path / "results-acmg-pick"
    env = dict(os.environ, PATH=f"{stub}:{os.environ['PATH']}")
    submit = HERE.parent / "submit_acmg_carriers.sh"
    res = subprocess.run(["bash", str(submit), "--batch-root", str(batch), "-P", "acc_x",
                          "--mode", "pick", "--out-root", str(out_root),
                          "--pick-order", "rank,mane_select"],
                         capture_output=True, text=True, env=env, cwd=tmp_path)
    assert res.returncode == 0, res.stderr
    for chrom in ("chr21", "chr22"):
        v = list(csv.DictReader(
            (out_root / f"results-carriers-pick-{chrom}_v6" / "acmg_variants.tsv").open(),
            delimiter="\t"))
        # rank first -> stop_gained ENST1 wins over the canonical missense
        assert [(r["transcript_id"], r["pick_tier"]) for r in v] == [("ENST1", "rank")]
        c = (out_root / f"results-carriers-pick-{chrom}_v6" / "acmg_carriers.tsv").read_text()
        assert "\themi\n" in c
    assert sorted(str(p) for p in batch.rglob("*")) == before

    # pick mode refuses a missing or nested --out-root
    for extra in ([], ["--out-root", str(batch / "x")]):
        bad = subprocess.run(["bash", str(submit), "--batch-root", str(batch), "-P", "acc_x",
                              "--mode", "pick", *extra], capture_output=True, text=True,
                             env=env, cwd=tmp_path)
        assert bad.returncode != 0
    assert sorted(str(p) for p in batch.rglob("*")) == before


def test_compare_mane_vs_pick(tmp_path):
    """Run both scripts on one synthetic chunk and compare: a gene with no MANE
    transcript gains a P/LP call under pick; a hemizygous carrier appears."""
    rec = [
        ("X", 100, [blk("stop_gained", "GNOMANE", "GN", "ENST1", "P", "PVS1", canonical="YES",
                        impact="HIGH")], ("1", "0/1", "0/0")),
        ("X", 200, [blk("missense_variant", "GM", "GM", "ENST2", "LP", "PM1", mane="NM_2")],
         ("0/1", "0/0", "0/0")),
    ]
    vcf = write_vcf(tmp_path, rec)
    mane_root, pick_root = tmp_path / "mane", tmp_path / "pick"
    md, pd = mane_root / "results-carriers-chrX_v6", pick_root / "results-carriers-pick-chrX_v6"
    md.mkdir(parents=True)
    pd.mkdir(parents=True)
    for script, d in (("build_acmg_carriers.py", md), ("build_acmg_carriers_pick.py", pd)):
        subprocess.run([sys.executable, str(HERE.parent / script), "--vcf", str(vcf),
                        "--out-variants", str(d / "acmg_variants.tsv"),
                        "--out-carriers", str(d / "acmg_carriers.tsv")],
                       check=True, capture_output=True)
    out = pick_root / "comparison"
    subprocess.run([sys.executable, str(HERE.parent / "compare_acmg_pick_vs_mane.py"),
                    "--mane-root", str(mane_root), "--pick-root", str(pick_root),
                    "--out", str(out)], check=True, capture_output=True)
    overall = {r["metric"]: r for r in csv.DictReader((out / "overall.tsv").open(), delimiter="\t")}
    assert overall["plp_unique_variants"]["mane_only"] == "1"
    assert overall["plp_unique_variants"]["pick"] == "2"
    assert overall["unique_carriers"]["pick"] == "2"
    assert overall["unique_hemizygous_carriers"]["pick"] == "1"
    genes = {r["gene"]: r for r in csv.DictReader((out / "per_gene.tsv").open(), delimiter="\t")}
    assert genes["GNOMANE"]["plp_variants_delta"] == "1"
    assert genes["GNOMANE"]["hemi_carriers_pick"] == "1"
    only = list(csv.DictReader((out / "plp_only_in_one.tsv").open(), delimiter="\t"))
    assert [(r["gene"], r["found_in"]) for r in only] == [("GNOMANE", "pick_only")]
    tiers = {r["pick_tier"]: r["rows"] for r in
             csv.DictReader((out / "pick_tier_counts.tsv").open(), delimiter="\t")}
    assert tiers == {"only_candidate": "2"}


def test_submit_only_exclude_and_limits(tmp_path):
    """--only chr1 must not also take chr10; --exclude drops a chromosome;
    -W/-n reach bsub (overriding the task file's #BSUB defaults)."""
    import os
    batch = tmp_path / "batch"
    for c in ("chr1", "chr10", "chr2"):
        (batch / f"results-fastvep-{c}_v6").mkdir(parents=True)
    stub = tmp_path / "bin"
    stub.mkdir()
    log = tmp_path / "bsub_args"
    (stub / "bsub").write_text(f'#!/usr/bin/env bash\nprintf "%s\\n" "$@" > {log}\ncat >/dev/null\n')
    (stub / "bsub").chmod(0o755)
    env = dict(os.environ, PATH=f"{stub}:{os.environ['PATH']}")
    submit = HERE.parent / "submit_acmg_carriers.sh"

    def sub(*extra):
        res = subprocess.run(["bash", str(submit), "--batch-root", str(batch), "-P", "acc_x",
                              "--mode", "pick", "--out-root", str(tmp_path / "out"), *extra],
                             capture_output=True, text=True, env=env, cwd=tmp_path)
        assert res.returncode == 0, res.stderr
        manifest = max((tmp_path / "out").glob(".acmg_carriers_manifest.*"),
                       key=lambda p: p.stat().st_mtime_ns)
        chrs = [Path(l).name for l in manifest.read_text().split()]
        manifest.unlink()
        return chrs, log.read_text().split("\n")

    chrs, args = sub("--only", "chr1", "-W", "24:00", "-n", "12")
    assert chrs == ["results-fastvep-chr1_v6"]
    assert args[args.index("-W") + 1] == "24:00" and args[args.index("-n") + 1] == "12"
    assert "AC_JOBS=12" in " ".join(args)

    chrs, args = sub("--exclude", "chr1")
    assert sorted(chrs) == ["results-fastvep-chr10_v6", "results-fastvep-chr2_v6"]
    assert "-W" not in args and "-n" not in args

    bad = subprocess.run(["bash", str(submit), "--batch-root", str(batch), "-P", "acc_x",
                          "--mode", "pick", "--out-root", str(tmp_path / "out"),
                          "--exclude", "chr99"], capture_output=True, text=True, env=env,
                         cwd=tmp_path)
    assert bad.returncode == 3
