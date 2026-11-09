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

    # pick mode refuses a missing --out-root, the batch root itself, and
    # anything inside a fastVEP chunk dir; nothing is created when it refuses.
    for extra in ([], ["--out-root", str(batch)],
                  ["--out-root", str(batch / "results-fastvep-chr21_v6" / "x")]):
        bad = subprocess.run(["bash", str(submit), "--batch-root", str(batch), "-P", "acc_x",
                              "--mode", "pick", *extra], capture_output=True, text=True,
                             env=env, cwd=tmp_path)
        assert bad.returncode != 0, extra
    assert sorted(str(p) for p in batch.rglob("*")) == before

    # ...but a separate subfolder of the batch root is allowed, and the run
    # writes only inside it.
    sub = batch / "results-acmg-pick"
    ok = subprocess.run(["bash", str(submit), "--batch-root", str(batch), "-P", "acc_x",
                         "--mode", "pick", "--out-root", str(sub)], capture_output=True,
                        text=True, env=env, cwd=tmp_path)
    assert ok.returncode == 0, ok.stderr
    assert (sub / "results-carriers-pick-chr21_v6" / "acmg_variants.tsv").exists()
    after = [p for p in sorted(str(p) for p in batch.rglob("*")) if not p.startswith(str(sub))]
    assert after == before


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
    # refuses to write into a MANE result dir
    bad = subprocess.run([sys.executable, str(HERE.parent / "compare_acmg_pick_vs_mane.py"),
                          "--mane-root", str(mane_root), "--pick-root", str(pick_root),
                          "--out", str(md / "cmp")], capture_output=True)
    assert bad.returncode != 0 and not (md / "cmp").exists()
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


FV_HEADERS = [
    '##INFO=<ID=FV_CLINVAR,Number=.,Type=String,Description="fastVEP ClinVar annotations. '
    'Format: ALLELE|SIGNIFICANCE|REVIEW_STATUS|PHENOTYPES|VARIANT_CLASS|SO_ACCESSION">',
    '##INFO=<ID=FV_GNOMAD,Number=.,Type=String,Description="fastVEP gnomAD annotations. Format: '
    'ALLELE|ALL_AF|ALL_AC|ALL_AN|ALL_HC|AFR_AF|AMR_AF|ASJ_AF|EAS_AF|FIN_AF|MID_AF|NFE_AF|OTH_AF|'
    'REMAINING_AF|SAS_AF">',
    '##INFO=<ID=FV_REVEL,Number=.,Type=String,Description="fastVEP REVEL annotations. '
    'Format: ALLELE|SCORE">',
    '##INFO=<ID=FV_GNOMAD_GENE,Number=.,Type=String,Description="fastVEP gnomAD gene constraint '
    'annotations. Format: SYMBOL|PLI|LOEUF|MIS_Z|SYN_Z">',
]
EXTRA_FIELDS = FIELDS + "|Protein_position|Amino_acids|CDS_position|FLAGS|Existing_variation"


def test_functional_and_fv_columns(tmp_path):
    """New columns: per-transcript protein fields + is_last_exon from the picked
    block; ClinVar/gnomAD/REVEL/constraint from FV_* INFO, matched by allele
    (VCF alt for SNVs, VEP-style '-' for a deletion) and by gene symbol."""
    def blk2(cons, sym, gene, feat, exon, extra, **kw):
        b = blk(cons, sym, gene, feat, **kw).split("|")
        b[7] = exon                       # EXON
        return "|".join(b + list(extra))

    header = ("##fileformat=VCFv4.2\n"
              f'##INFO=<ID=CSQ,Number=.,Type=String,Description="fastVEP. Format: {EXTRA_FIELDS}">\n'
              + "\n".join(FV_HEADERS) + "\n"
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    # 1) the real chr21 record's INFO (LIPI), SNV, alt allele "A"
    snv_csq = blk2("missense_variant", "LIPI", "G_LIPI", "ENST1", "3/10",
                   ("120", "R/S", "358", "", "rs123"), acmg="VUS", mane="NM_1", allele="A")
    snv_info = ("AF=1.2e-05;"
                "FV_CLINVAR=A|Uncertain_significance|criteria_provided%2C_single_submitter|"
                "not_specified|single_nucleotide_variant|SO%3A0001483;"
                "FV_GNOMAD=A|1.5e-6|1|627454|0|0.0000565|0|0|0|0|0|0||0|0;"
                "FV_PHYLOP=A|-0.438;FV_REVEL=A|0.17;FV_OMIM=LIPI|609252|;"
                "FV_GNOMAD_GENE=LIPI|0.0|1.875|-0.31|-0.71,OTHER|1.0|0.1|0|0;"
                f"CSQ={snv_csq}")
    # 2) a deletion in the last exon of a gene with an incomplete CDS; FV
    #    fields keyed by the VEP-style allele "-"
    del_csq = blk2("frameshift_variant", "GD", "G_D", "ENST2", "12/12",
                   ("455", "K/X", "1364", "cds_end_NF", ""), acmg="LP", crit="PVS1",
                   canonical="YES", impact="HIGH", allele="-")
    del_info = ("FV_GNOMAD=-|0.001|5|5000|1|0.002|0.0001|0.05|||||||0.003;"
                "FV_GNOMAD_GENE=GD|0.99|0.2|1|1;"
                f"CSQ={del_csq}")
    # 3) FV fields present but for a different allele -> counted as unmatched
    odd_info = "FV_REVEL=G|0.9;CSQ=" + blk2("missense_variant", "GO", "G_O", "ENST3", "",
                                             ("", "", "", "", ""), allele="T")
    body = "\n".join([f"21\t14109002\t.\tC\tA\t.\t.\t{snv_info}",
                      f"1\t500\t.\tCA\tC\t.\t.\t{del_info}",
                      f"1\t900\t.\tC\tT\t.\t.\t{odd_info}"]) + "\n"
    vcf = tmp_path / "fv.vcf"
    vcf.write_text(header + body, encoding="utf-8")
    ov, oc = tmp_path / "v.tsv", tmp_path / "c.tsv"
    res = subprocess.run([sys.executable, str(SCRIPT), "--vcf", str(vcf), "--out-variants",
                          str(ov), "--out-carriers", str(oc)],
                         check=True, capture_output=True, text=True)
    rows = list(csv.DictReader(ov.open(), delimiter="\t"))
    assert list(rows[0].keys())[:26] == list(P.VARIANT_COLUMNS[:26])  # old columns unchanged
    by = {r["gene"]: r for r in rows}

    r = by["LIPI"]
    assert (r["protein_position"], r["amino_acids"], r["cds_position"], r["existing_variation"]) \
        == ("120", "R/S", "358", "rs123")
    assert r["is_last_exon"] == "0"
    assert r["clinvar_significance"] == "Uncertain_significance"
    assert r["clinvar_review_status"] == "criteria_provided,_single_submitter"
    assert (r["gnomad_af"], r["gnomad_grpmax_af"], r["gnomad_grpmax_pop"], r["gnomad_hom_count"]) \
        == ("1.5e-6", "0.0000565", "AFR", "0")
    assert r["revel"] == "0.17"
    assert (r["gene_pli"], r["gene_loeuf"]) == ("0.0", "1.875")   # LIPI's, not OTHER's

    r = by["GD"]
    assert r["is_lof_consequence"] == "1" and r["is_last_exon"] == "1"
    assert r["flags"] == "cds_end_NF"
    # ASJ (0.05) is excluded from grpmax; SAS 0.003 is the highest kept group
    assert (r["gnomad_af"], r["gnomad_grpmax_af"], r["gnomad_grpmax_pop"]) == ("0.001", "0.003", "SAS")
    assert r["gnomad_hom_count"] == "1" and r["clinvar_significance"] == ""
    assert (r["gene_pli"], r["gene_loeuf"]) == ("0.99", "0.2")

    r = by["GO"]
    assert r["revel"] == "" and r["is_last_exon"] == ""
    assert "had no entry for the alt allele: 1" in res.stderr


def test_vep_allele_and_last_exon():
    assert P.vep_allele("C", "A") == "A"
    assert P.vep_allele("CA", "C") == "-"
    assert P.vep_allele("C", "CTT") == "TT"
    assert P.vep_allele("CAT", "GT") == "GT"           # no shared anchor: unchanged
    assert P.is_last_exon("12/12") == "1" and P.is_last_exon("11-12/12") == "1"
    assert P.is_last_exon("3/10") == "0" and P.is_last_exon("") == ""


def test_extract_gene_carriers(tmp_path):
    """Carriers of listed genes only, across chromosomes, annotation joined,
    compound-het count per (person, gene)."""
    root = tmp_path / "pick"
    recs = {
        "chr1_v6": [("1", 100, [blk("stop_gained", "MUTYH", "GM", "ENSTM", "P", "PVS1",
                                     mane="NM_M", impact="HIGH")], ("0/1", "0/1", "0/0")),
                    ("1", 200, [blk("frameshift_variant", "MUTYH", "GM", "ENSTM", "LP", "PVS1",
                                     mane="NM_M", impact="HIGH")], ("0/1", "0/0", "0/0")),
                    ("1", 300, [blk("stop_gained", "OTHER", "GO", "ENSTO", "P", "PVS1",
                                     mane="NM_O", impact="HIGH")], ("1/1", "0/0", "0/0"))],
        "chr17_v6": [("17", 400, [blk("stop_gained", "TP53", "GT", "ENSTT", "P", "PVS1",
                                       mane="NM_T", impact="HIGH")], ("0/0", "0/0", "1/1"))],
    }
    for chrom, rr in recs.items():
        d = root / f"results-carriers-pick-{chrom}"
        d.mkdir(parents=True)
        sub = tmp_path / chrom
        sub.mkdir()
        vcf = write_vcf(sub, rr)
        subprocess.run([sys.executable, str(SCRIPT), "--vcf", str(vcf),
                        "--out-variants", str(d / "acmg_variants.tsv"),
                        "--out-carriers", str(d / "acmg_carriers.tsv")],
                       check=True, capture_output=True)
    genes = tmp_path / "genes.txt"
    genes.write_text("# test list\nMUTYH\nTP53, BRCA1\nNOTAGENE\n")
    out, summ = tmp_path / "carriers.tsv", tmp_path / "summary.tsv"
    res = subprocess.run([sys.executable, str(HERE.parent / "extract_gene_carriers.py"),
                          "--pick-root", str(root), "--genes", str(genes), "--out", str(out),
                          "--out-summary", str(summ)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    got = sorted((r["gene"], r["pos"], r["person_id"], r["n_plp_variants_in_gene"]) for r in rows)
    assert got == [("MUTYH", "100", "S1", "2"), ("MUTYH", "100", "S2", "1"),
                   ("MUTYH", "200", "S1", "2"), ("TP53", "400", "S3", "1")]
    assert rows[0]["acmg_label"] in ("Pathogenic", "Likely_pathogenic")
    assert "transcript_id" in rows[0] and "gnomad_af" in rows[0]
    s = {r["gene"]: r for r in csv.DictReader(summ.open(), delimiter="\t")}
    assert (s["MUTYH"]["plp_variants_with_carriers"], s["MUTYH"]["carriers"],
            s["MUTYH"]["carriers_2plus_plp_variants"]) == ("2", "2", "1")
    assert s["TP53"]["hom_carriers"] == "1" and s["BRCA1"]["carriers"] == "0"
    assert "NOTAGENE" in res.stderr and "BRCA1" in res.stderr
