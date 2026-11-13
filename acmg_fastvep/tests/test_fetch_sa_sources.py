"""fetch_sa_sources.sh against local fixtures, with curl stubbed.

The stub behaves like `curl -f` against a real server in the one way that
matters here: resuming (`-C -`) a file that is already complete gets HTTP 416,
which `-f` turns into exit 22. Two bugs this file pins down:

* `gzip -dc big.vcf.gz | head -300 | grep -q ...` under `set -o pipefail` fails
  on every file big enough that gzip is still writing when head stops reading
  (SIGPIPE), so a perfectly good ClinVar VCF was rejected on the first run.
* a re-run over a file that finished downloading but was never marked verified
  resumed it, got 416 / exit 22, and could never get past that file again.

Needs tabix + bgzip (htslib) and unzip; skipped without them.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import random
import shutil
import subprocess
import zipfile
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "fetch_sa_sources.sh"

pytestmark = pytest.mark.skipif(
    not all(shutil.which(t) for t in ("tabix", "bgzip", "unzip")),
    reason="needs tabix, bgzip and unzip")

STUB_CURL = r"""#!/usr/bin/env bash
out=""; resume=0; wfmt=""; url=""
while [ $# -gt 0 ]; do
  case "$1" in
    -o) out="$2"; shift 2 ;;
    -C) resume=1; shift 2 ;;
    -w) wfmt="$2"; shift 2 ;;
    --retry|--retry-delay) shift 2 ;;
    http*) url="$1"; shift ;;
    *) shift ;;
  esac
done
fix="$FIXTURES/$(basename "$url")"
say() { [ -n "$wfmt" ] && printf '%s' "$1"; }
if [ ! -f "$fix" ]; then echo "GET $url 404" >> "$CURL_LOG"; say 404; exit 22; fi
size() { wc -c < "$1" | tr -d ' '; }
if [ "$resume" = 1 ] && [ -f "$out" ]; then
  have=$(size "$out"); want=$(size "$fix")
  if [ "$have" -ge "$want" ]; then echo "GET $url 416" >> "$CURL_LOG"; say 416; exit 22; fi
  tail -c +"$((have + 1))" "$fix" >> "$out"
  echo "GET $url 206" >> "$CURL_LOG"; say 206; exit 0
fi
cp "$fix" "$out"
echo "GET $url 200" >> "$CURL_LOG"; say 200; exit 0
"""


def _gz(path: Path, text: str):
    with gzip.open(path, "wt") as fh:
        fh.write(text)


@pytest.fixture
def env(tmp_path):
    fix = tmp_path / "fixtures"
    fix.mkdir()
    # ClinVar: big enough that gzip is still writing when the header check stops reading.
    hdr = ("##fileformat=VCFv4.1\n##fileDate=2026-11-09\n"
           '##INFO=<ID=AF_EXAC,Number=1,Type=Float,Description="x">\n'
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    rnd = random.Random(1)
    body = "".join(f"1\t{i}\t{rnd.randrange(10**9)}\tA\tG\t.\t.\tAF_EXAC=0.{rnd.randrange(10**6)}\n"
                   for i in range(1, 300_001))
    for name in ("clinvar_20261109.vcf.gz", "clinvar.vcf.gz"):
        _gz(fix / name, hdr + body)
        md5 = hashlib.md5((fix / name).read_bytes()).hexdigest()
        (fix / (name + ".md5")).write_text(f"{md5}  {name}\n")
    _gz(fix / "variant_summary.txt.gz", "#AlleleID\tType\n1\tsnv\n")
    # gnomAD chr21, really bgzipped and indexed
    g = fix / "gnomad.exomes.v4.1.sites.chr21.vcf"
    g.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 + "".join(f"chr21\t{1000 + i}\t.\tA\tG\t.\tPASS\tAC=1\n" for i in range(20000)))
    subprocess.run(["bgzip", "-f", str(g)], check=True)
    os.rename(str(g) + ".gz", str(g) + ".bgz")
    subprocess.run(["tabix", "-p", "vcf", str(g) + ".bgz"], check=True)
    (fix / "gnomad.v4.1.constraint_metrics.tsv").write_text(
        "gene\tpLI\n" + "".join(f"G{i}\t0.5\n" for i in range(10001)))
    with zipfile.ZipFile(fix / "revel-v1.3_all_chromosomes.zip", "w") as z:
        z.writestr("revel_with_transcript_ids",
                   "chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid\n"
                   "1,10,10,A,G,M,V,0.5,ENST1\n21,20,20,C,T,R,W,0.9,ENST2\n")
    (fix / "download").write_text('"GENE SYMBOL","DISEASE LABEL"\n'
                                  + "".join(f'"G{i}","d"\n' for i in range(1001)))
    _gz(fix / "rmsk.txt.gz", "0\t1\t2\n")

    stub = tmp_path / "bin"
    stub.mkdir()
    (stub / "curl").write_text(STUB_CURL)
    (stub / "curl").chmod(0o755)
    home = tmp_path / "fv_home"
    home.mkdir()
    e = {k: v for k, v in os.environ.items() if not k.startswith("FV_")}
    e.update({"PATH": f"{stub}:{e['PATH']}", "FIXTURES": str(fix), "FV_HOME": str(home),
              "CURL_LOG": str(tmp_path / "curl.log")})
    return {"env": e, "src": home / "sa_src", "log": tmp_path / "curl.log", "fix": fix}


def _run(env, *args):
    return subprocess.run(["bash", str(SCRIPT), "--chroms", "21", *args], env=env["env"],
                          capture_output=True, text=True)


def _gets(env, what):
    return [ln for ln in env["log"].read_text().splitlines() if what in ln] if env["log"].exists() else []


def test_first_run_accepts_a_large_clinvar(env):
    """The AF_EXAC header check must pass on a ClinVar VCF larger than a pipe buffer."""
    r = _run(env, "--clinvar-date", "20261109")
    assert r.returncode == 0, r.stderr
    src = env["src"]
    assert (src / "clinvar_20261109.vcf.gz.md5.ok").exists()
    assert (src / "CLINVAR_DATE").read_text().strip() == "20261109"
    assert (src / "VARIANT_SUMMARY_DATE").read_text().strip().isdigit()
    assert (src / "revel" / "revel_chr21.csv").exists()
    rows = (src / "SOURCES.tsv").read_text().splitlines()
    assert rows[0].split("\t")[:5] == ["file", "url", "date", "size", "sha256"]
    names = {ln.split("\t")[0] for ln in rows[1:]}
    assert {"clinvar_20261109.vcf.gz", "variant_summary.txt.gz", "gnomad.exomes.v4.1.sites.chr21.vcf.bgz",
            "gnomad.v4.1.constraint_metrics.tsv", "revel-v1.3_all_chromosomes.zip",
            "rmsk.txt.gz"} <= names
    assert any("MOVING file" in ln for ln in rows if ln.startswith("variant_summary"))


def test_latest_clinvar_is_named_by_its_filedate(env):
    r = _run(env)
    assert r.returncode == 0, r.stderr
    assert (env["src"] / "clinvar_20261109.vcf.gz").exists()
    assert not (env["src"] / "clinvar.latest.vcf.gz").exists()


def test_rerun_over_complete_but_unmarked_files(env):
    """Complete downloads whose .ok marker was never written (the run died after
    the download) must be accepted on the re-run, not resumed into a 416."""
    assert _run(env, "--clinvar-date", "20261109").returncode == 0
    src = env["src"]
    for marker in ("clinvar_20261109.vcf.gz.md5.ok", "gnomad.exomes.v4.1.sites.chr21.vcf.bgz.ok",
                   "revel/.ok"):
        (src / marker).unlink()
    env["log"].unlink()
    r = _run(env, "--clinvar-date", "20261109")
    assert r.returncode == 0, r.stderr
    for big in ("clinvar_20261109.vcf.gz", "gnomad.exomes.v4.1.sites.chr21.vcf.bgz",
                "revel-v1.3_all_chromosomes.zip"):
        assert not [g for g in _gets(env, big) if not g.split()[1].endswith((".md5", ".tbi"))], \
            f"{big} was downloaded again"
    assert "already complete" in r.stderr


def test_416_on_resume_is_success_then_verified(env):
    """If the integrity check is not what said 'complete', a resume that gets 416
    is still not an error: the file is checked and accepted."""
    assert _run(env, "--clinvar-date", "20261109").returncode == 0
    g = env["src"] / "gnomad.exomes.v4.1.sites.chr21.vcf.bgz"
    (env["src"] / (g.name + ".ok")).unlink()
    (env["src"] / (g.name + ".tbi")).unlink()     # gnomad_ok fails -> fetch path
    env["log"].unlink()
    r = _run(env, "--clinvar-date", "20261109")
    assert r.returncode == 0, r.stderr
    assert any(g.name + " 416" in ln for ln in _gets(env, g.name))
    assert "HTTP 416 on resume" in r.stderr


def test_truncated_gnomad_is_resumed(env):
    assert _run(env, "--clinvar-date", "20261109").returncode == 0
    g = env["src"] / "gnomad.exomes.v4.1.sites.chr21.vcf.bgz"
    (env["src"] / (g.name + ".ok")).unlink()
    data = g.read_bytes()
    g.write_bytes(data[: len(data) // 2])
    env["log"].unlink()
    r = _run(env, "--clinvar-date", "20261109")
    assert r.returncode == 0, r.stderr
    assert any(" 206" in ln for ln in _gets(env, g.name + " "))
    assert g.read_bytes() == data


def test_bad_md5_still_fails(env):
    (env["fix"] / "clinvar_20261109.vcf.gz.md5").write_text("0" * 32 + "  x\n")
    r = _run(env, "--clinvar-date", "20261109")
    assert r.returncode != 0 and "md5" in r.stderr
