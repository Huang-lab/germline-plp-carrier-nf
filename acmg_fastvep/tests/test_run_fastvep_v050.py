"""run_fastvep.sh / run_fastvep_batch.sh / run_fastvep.lsf against the v0.5.0 contract.

fastvep is a stub that records its argv and can print one of the stderr lines
fastVEP v0.5.0 uses to say the reference stack is not what it should be (it
still exits 0 then). Only bash and python3 are needed.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ACMG = HERE.parent
RUNNER = ACMG / "run_fastvep.sh"
BATCH = ACMG / "run_fastvep_batch.sh"
LSF = ACMG / "run_fastvep.lsf"

# The four fastVEP v0.5.0 messages run_fastvep.sh gates on (annotate.rs / lib.rs at the tag).
GATED = [
    "warning: --acmg is running on an incomplete reference stack: REVEL (revel).",
    "/x/sa_db/revel_chr1.osa2 was not loaded: it holds no record",
    "warning: --acmg: the gnomAD databases that are loaded hold nothing of chromosomes 22, "
    "which the gene model has transcripts on.",
    "warning: 3 allele lookup(s) were answered by more than one database under the key `gnomad`",
]


def _fake(tmp_path: Path, *, stderr_line: str = "", tag: str = "a") -> Path:
    p = tmp_path / f"fake_fastvep_{tag}"
    p.write_text(
        "#!/usr/bin/env bash\n"
        f"# {tag}\n"
        '[ "$1" = --version ] && { echo "fastvep 0.5.0"; exit 0; }\n'
        'echo "$*" >> "$ARGV_LOG"\n'
        'printf "##fileformat=VCFv4.2\\n"\n'
        'printf "##INFO=<ID=CSQ,Number=.,Type=String,Description=\\"fastVEP. '
        'Format: Allele|SYMBOL|ACMG|ACMG_CRITERIA|PICK|ACMG_DISORDERS\\">\\n"\n'
        'printf "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n"\n'
        'printf "17\\t100\\t.\\tA\\tG\\t.\\tPASS\\tCSQ=G|BRCA1|LP|PVS1|1|\\n"\n'
        + (f"echo {stderr_line!r} >&2\n" if stderr_line else "")
        + "exit 0\n"
    )
    p.chmod(0o755)
    return p


def _setup(tmp_path: Path, fake: Path, *, sealed_sa: bool = True):
    indir = tmp_path / "norm_qc"
    indir.mkdir(exist_ok=True)
    vcf = indir / "chunk1.norm.vcf.gz"
    with gzip.open(vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 "17\t100\t.\tA\tG\t.\tPASS\t.\n")
    sa = tmp_path / "sa_db"
    sa.mkdir(exist_ok=True)
    if sealed_sa:
        (sa / "SA_DB_SHA256.txt").write_text("abc  clinvar.osa2\n")
    for name, body in (("g.gff3", "##gff-version 3\n"), ("g.fa", ">17\nA\n"), ("g.cache", "c\n")):
        (tmp_path / name).write_text(body)
    env = {k: v for k, v in os.environ.items() if not k.startswith("FV_")}
    env.update({
        "FV_BIN": str(fake), "FV_BIN_SHA256": hashlib.sha256(fake.read_bytes()).hexdigest(),
        "FV_GFF3": str(tmp_path / "g.gff3"), "FV_FASTA": str(tmp_path / "g.fa"),
        "FV_CACHE": str(tmp_path / "g.cache"), "FV_SA": str(sa),
        "ARGV_LOG": str(tmp_path / "argv.log"),
    })
    return vcf, env


def _run(args, env):
    return subprocess.run(["bash", *map(str, args)], env=env, capture_output=True, text=True)


def test_acmg_command_line(tmp_path):
    fake = _fake(tmp_path)
    vcf, env = _setup(tmp_path, fake)
    out = tmp_path / "out"
    r = _run([RUNNER, "-i", vcf, "-o", out, "--acmg"], env)
    assert r.returncode == 0, r.stderr
    argv = (tmp_path / "argv.log").read_text().split()
    for flag in ("--acmg", "--flag-pick-allele-gene", "--transcript-cache", "--fasta",
                 "--no-progress", "--sa-dir"):
        assert flag in argv, flag
    assert argv[argv.index("--transcript-cache") + 1] == env["FV_CACHE"]
    assert "--gff3" not in argv and "--pick" not in argv
    prov = (out / "chunk1.fastvep.provenance").read_text()
    assert f"fastvep_sha256   {env['FV_BIN_SHA256']}" in prov
    assert "sa_db_sha256     " + hashlib.sha256(
        (tmp_path / "sa_db" / "SA_DB_SHA256.txt").read_bytes()).hexdigest() in prov
    assert (out / "chunk1.acmg_plp.tsv").exists()
    assert (out / "chunk1.fastvep.stderr.log").exists()


def test_pick_and_acmg_are_exclusive(tmp_path):
    fake = _fake(tmp_path)
    vcf, env = _setup(tmp_path, fake)
    r = _run([RUNNER, "-i", vcf, "-o", tmp_path / "out", "--acmg", "--pick"], env)
    assert r.returncode == 3 and "--pick cannot be combined with --acmg" in r.stderr


@pytest.mark.parametrize("line", GATED)
def test_gated_stderr_fails_the_chunk(tmp_path, line):
    fake = _fake(tmp_path, stderr_line=line)
    vcf, env = _setup(tmp_path, fake)
    out = tmp_path / "out"
    r = _run([RUNNER, "-i", vcf, "-o", out, "--acmg"], env)
    assert r.returncode == 5, r.stderr
    assert not (out / "chunk1.fastvep.vcf.gz").exists()
    assert not (out / "chunk1.acmg_plp.tsv").exists()
    assert not (out / "chunk1.fastvep.provenance").exists()
    assert line in (out / "chunk1.fastvep.stderr.log").read_text()
    assert not list(out.glob(".*partial"))


def test_unpinned_binary_and_missing_env_refused(tmp_path):
    fake = _fake(tmp_path)
    vcf, env = _setup(tmp_path, fake)
    r = _run([RUNNER, "-i", vcf, "-o", tmp_path / "out"], {**env, "FV_BIN_SHA256": "0" * 64})
    assert r.returncode == 3 and "not the pinned binary" in r.stderr
    bare = {k: v for k, v in env.items() if k not in ("FV_CACHE", "FV_FASTA")}
    r = _run([RUNNER, "-i", vcf, "-o", tmp_path / "out"], bare)
    assert r.returncode == 3 and "FV_FASTA" in r.stderr and "FV_CACHE" in r.stderr
    sa_unsealed = tmp_path / "sa2"
    sa_unsealed.mkdir()
    r = _run([RUNNER, "-i", vcf, "-o", tmp_path / "out", "--acmg", "--sa-dir", sa_unsealed], env)
    assert r.returncode == 3 and "check_sa_db.sh" in r.stderr


def test_batch_skips_only_matching_provenance(tmp_path):
    fake = _fake(tmp_path, tag="a")
    _, env = _setup(tmp_path, fake)
    out = tmp_path / "out"
    args = [BATCH, "--in-dir", tmp_path / "norm_qc", "-o", out]
    r = _run(args, env)
    assert r.returncode == 0, r.stderr
    runprov = (out / "fastvep_run_provenance.txt").read_text()
    for key in ("fastvep_sha256", "transcript_cache_sha256", "sa_db_sha256", "pipeline_rev", "args"):
        assert key in runprov, key
    assert "\tok" in (out / "batch_manifest.tsv").read_text()

    r = _run(args, env)
    assert r.returncode == 0 and "skipped-already-done" in (out / "batch_manifest.tsv").read_text()

    # A different (but correctly pinned) binary must not skip over the old output.
    other = _fake(tmp_path, tag="b")
    env_b = {**env, "FV_BIN": str(other),
             "FV_BIN_SHA256": hashlib.sha256(other.read_bytes()).hexdigest()}
    r = _run(args, env_b)
    assert r.returncode == 6, r.stderr
    assert "chunk1: made by fastvep sha256" in r.stderr

    r = _run([*args, "--force"], env_b)
    assert r.returncode == 0, r.stderr
    assert env_b["FV_BIN_SHA256"] in (out / "chunk1.fastvep.provenance").read_text()

    # Output with no provenance at all (made by the pre-v0.5.0 scripts) is refused too.
    (out / "chunk1.fastvep.provenance").unlink()
    r = _run(args, env_b)
    assert r.returncode == 6 and "no provenance file" in r.stderr


def test_batch_acmg_needs_a_sealed_sa_db(tmp_path):
    fake = _fake(tmp_path)
    _, env = _setup(tmp_path, fake, sealed_sa=False)
    r = _run([BATCH, "--in-dir", tmp_path / "norm_qc", "-o", tmp_path / "out", "--acmg"], env)
    assert r.returncode == 3 and "check_sa_db.sh --quick failed" in r.stderr
    assert not (tmp_path / "argv.log").exists()


def test_lsf_requires_env(tmp_path):
    env = {k: v for k, v in os.environ.items() if not k.startswith("FV_")}
    env.update({"FV_IN": str(tmp_path), "FV_OUT": str(tmp_path / "o")})
    r = _run([LSF], env)
    assert r.returncode != 0 and "FV_ENV" in r.stderr


@pytest.mark.parametrize("script", [RUNNER, BATCH, LSF])
def test_scripts_are_syntactically_valid(script):
    assert subprocess.run(["bash", "-n", str(script)]).returncode == 0
