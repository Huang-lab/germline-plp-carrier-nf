"""Functional tests for build_sa_dbs.sh and check_sa_db.sh (fastVEP v0.5.0 stack).

fastvep and the fastVEP helper scripts are stubbed: the contract under test is
which databases get built from which inputs with which flags, that nothing
reaches sa_db unless sa-build succeeded, that `.done` gates the skip, and that
check_sa_db.sh accepts exactly the expected file set and nothing else.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
BUILD = HERE.parent / "build_sa_dbs.sh"
CHECK = HERE.parent / "check_sa_db.sh"
LSF = HERE.parent / "build_sa_dbs.lsf"

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]

# Stub fastvep: `sa-build` writes <-o>.<ext> as v0.5.0's --format auto would
# (gene sources .oga, custom_bed .osi, everything else .osa2), logs its argv,
# and fails for $FAIL_SOURCE.
STUB_FASTVEP = r"""#!/usr/bin/env bash
[ "$1" = --version ] && { echo "fastvep 0.5.0"; exit 0; }
echo "$*" >> "$STUB_LOG"
src=""; out=""
while [ $# -gt 0 ]; do
  case "$1" in --source) src="$2"; shift 2;; -o) out="$2"; shift 2;; *) shift;; esac
done
[ -n "${FAIL_SOURCE:-}" ] && [ "$src" = "$FAIL_SOURCE" ] && { echo "boom" >&2; exit 1; }
case "$src" in
  omim|gnomad_genes|clinvar_protein) ext=oga ;;
  custom_bed) ext=osi ;;
  *) ext=osa2 ;;
esac
echo "db $src" > "$out.$ext"
"""

STUB_EXTRACT = r"""import argparse, os
p = argparse.ArgumentParser()
p.add_argument("--chrom"); p.add_argument("--url"); p.add_argument("--out-dir")
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
open(os.path.join(a.out_dir, f"spliceai_chr{a.chrom}.vcf"), "w").write("x\n")
open(os.path.join(a.out_dir, f"phylop_chr{a.chrom}.tsv"), "w").write("x\n")
"""

STUB_CLINGEN = r"""import sys
src, dst = sys.argv[1], sys.argv[2]
with open(dst, "w") as f:
    for g in ("BRCA1", "BRCA2", "TP53", "PTEN", "MYH7"):
        cols = [""] * 13; cols[5] = "0"; cols[8] = g; cols[12] = "x (ClinGen Definitive/AD, MONDO:1)"
        f.write("\t".join(cols) + "\n")
"""

STUB_RMSK = r"""import sys
print("chr1\t0\t10\tAluY")
"""

# check_acmg_stack.py stub: writes a probe JSON whose BRCA1/PTEN/MYH7 genes carry `omim`.
STUB_STACK = r"""import argparse, json, os, sys
p = argparse.ArgumentParser()
for a in ("--gff3", "--fasta", "--sa-dir", "--fastvep", "--keep-json"): p.add_argument(a)
a = p.parse_args()
if os.environ.get("STACK_FAIL"):
    print("INCOMPLETE"); sys.exit(1)
ids = ["BRCA1_missense", "BRCA1_nonsense", "PTEN_missense", "MYH7_missense"]
omim = {} if os.environ.get("STACK_NO_OMIM") else {"omim": {}}
json.dump([{"id": i, "genes": {"G": dict(omim)}} for i in ids], open(a.keep_json, "w"))
print("Complete")
"""


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.fixture
def fv(tmp_path):
    home = tmp_path / "fv_home"
    src = home / "sa_src"
    (src / "revel").mkdir(parents=True)
    fvsrc = tmp_path / "fv_src"
    (fvsrc / "scripts").mkdir(parents=True)
    sasrc = fvsrc / "analysis/acmg_benchmark/scripts/sa_sources"
    sasrc.mkdir(parents=True)
    (fvsrc / "scripts/extract_gnomad_scores.py").write_text(STUB_EXTRACT)
    (fvsrc / "scripts/check_acmg_stack.py").write_text(STUB_STACK)
    (sasrc / "clingen_gdv_to_oga.py").write_text(STUB_CLINGEN)
    (sasrc / "repeatmasker_to_bed.py").write_text(STUB_RMSK)
    binp = tmp_path / "fastvep"
    binp.write_text(STUB_FASTVEP)
    binp.chmod(0o755)

    # sa_src as fetch_sa_sources.sh leaves it
    (src / "CLINVAR_DATE").write_text("20261109\n")
    for name in ("clinvar_20261109.vcf.gz", "clinvar_20261109.vcf.gz.md5.ok", "variant_summary.txt.gz",
                 "gnomad.v4.1.constraint_metrics.tsv", "clingen_gene_validity_20261110.csv",
                 "rmsk.txt.gz"):
        (src / name).write_text("data\n")
    for c in CHROMS:
        g = src / f"gnomad.exomes.v4.1.sites.chr{c}.vcf.bgz"
        g.write_text("data\n")
        Path(str(g) + ".tbi").write_text("idx\n")
        Path(str(g) + ".sha256").write_text("abc123\n")
        if c != "Y":
            (src / "revel" / f"revel_chr{c}.csv").write_text("h\nrow\n")
    gff3 = tmp_path / "g.gff3"; gff3.write_text("x\n")
    fasta = tmp_path / "g.fa"; fasta.write_text(">1\nA\n")
    env = {
        "PATH": "/usr/bin:/bin:/usr/local/bin:/opt/homebrew/bin",
        "FV_HOME": str(home), "FV_BIN": str(binp), "FV_BIN_SHA256": _sha(binp),
        "FV_SRC": str(fvsrc), "FV_SA": str(home / "sa_db"), "FV_GFF3": str(gff3),
        "FV_FASTA": str(fasta), "STUB_LOG": str(tmp_path / "stub.log"),
    }
    return {"env": env, "home": home, "db": home / "sa_db", "src": src, "log": tmp_path / "stub.log"}


def _run(script, fv, *args, extra=None):
    return subprocess.run(["bash", str(script), *args], env={**fv["env"], **(extra or {})},
                          capture_output=True, text=True)


def _calls(fv):
    return fv["log"].read_text().splitlines() if fv["log"].exists() else []


def test_single_builds_the_five_whole_genome_dbs(fv):
    r = _run(BUILD, fv, "--single")
    assert r.returncode == 0, r.stderr
    for f in ("clinvar.osa2", "clinvar_protein.oga", "gnomad_genes.oga", "omim.oga", "repeatmasker.osi"):
        assert (fv["db"] / f).stat().st_size > 0, f
        assert (fv["db"] / (f + ".done")).exists(), f
    calls = _calls(fv)
    assert all("--no-progress" in c and "--assembly GRCh38" in c for c in calls)
    assert not any("--allow-empty" in c for c in calls)
    by_src = {c.split("--source ")[1].split()[0]: c for c in calls}
    assert "clinvar_20261109.vcf.gz" in by_src["clinvar"]
    assert "variant_summary.txt.gz" in by_src["clinvar_protein"]     # not the VCF
    assert "derived/clingen_gdv.tsv" in by_src["omim"]               # ClinGen, not OMIM
    assert "--name repeatmasker" in by_src["custom_bed"]
    # nothing half-built is left in sa_db
    assert sorted(p.name for p in fv["db"].iterdir() if not p.name.endswith(".done")) == \
        ["clinvar.osa2", "clinvar_protein.oga", "gnomad_genes.oga", "omim.oga", "repeatmasker.osi"]


def test_use_omim_reads_genemap2(fv):
    (fv["src"] / "genemap2.txt").write_text("data\n")
    r = _run(BUILD, fv, "--single", "--use-omim")
    assert r.returncode == 0, r.stderr
    assert any("--source omim" in c and "genemap2.txt" in c for c in _calls(fv))


def test_chrom_builds_four_with_chunk_bits(fv):
    r = _run(BUILD, fv, "--chrom", "21")
    assert r.returncode == 0, r.stderr
    calls = {c.split("--source ")[1].split()[0]: c for c in _calls(fv)}
    assert set(calls) == {"gnomad", "spliceai", "phylop", "revel"}
    for s in ("gnomad", "spliceai", "phylop"):
        assert "--chunk-bits 16" in calls[s], calls[s]
    assert "--chunk-bits" not in calls["revel"]
    assert "gnomad.exomes.v4.1.sites.chr21.vcf.bgz" in calls["gnomad"]
    assert "spliceai_chr21.vcf" in calls["spliceai"] and "phylop_chr21.tsv" in calls["phylop"]
    for f in ("gnomad_chr21", "spliceai_chr21", "phylop_chr21", "revel_chr21"):
        assert (fv["db"] / f"{f}.osa2.done").exists()
    done = (fv["db"] / "gnomad_chr21.osa2.done").read_text()
    g = fv["src"] / "gnomad.exomes.v4.1.sites.chr21.vcf.bgz"
    assert f"source_sha256 {_sha(g)}" in done and f"source_file   {g}" in done
    # spliceai/phylop are keyed on the gnomAD file they were distilled from
    assert f"source_sha256 {_sha(g)}" in (fv["db"] / "spliceai_chr21.osa2.done").read_text()


def test_chrom_y_has_no_revel(fv):
    r = _run(BUILD, fv, "--chrom", "Y")
    assert r.returncode == 0, r.stderr
    assert not any("--source revel" in c for c in _calls(fv))
    assert "REVEL v1.3 has no chrY" in r.stderr


def test_skip_needs_output_and_done(fv):
    assert _run(BUILD, fv, "--chrom", "21").returncode == 0
    n = len(_calls(fv))
    assert _run(BUILD, fv, "--chrom", "21").returncode == 0
    assert len(_calls(fv)) == n, "re-run rebuilt something that had .done"
    (fv["db"] / "revel_chr21.osa2.done").unlink()
    r = _run(BUILD, fv, "--chrom", "21")
    assert r.returncode == 0
    new = _calls(fv)[n:]
    assert len(new) == 1 and "--source revel" in new[0]


def test_failed_build_publishes_nothing(fv):
    r = _run(BUILD, fv, "--chrom", "21", extra={"FAIL_SOURCE": "spliceai"})
    assert r.returncode != 0
    assert "spliceai" in r.stderr and "nothing published" in r.stderr
    assert not (fv["db"] / "spliceai_chr21.osa2").exists()
    assert not (fv["db"] / "spliceai_chr21.osa2.done").exists()
    assert (fv["db"] / "gnomad_chr21.osa2.done").exists()     # built before the failure


def test_missing_input_is_an_error(fv):
    (fv["src"] / "variant_summary.txt.gz").unlink()
    r = _run(BUILD, fv, "--single")
    assert r.returncode != 0 and "clinvar_protein" in r.stderr


def test_wrong_binary_hash_refused(fv):
    r = _run(BUILD, fv, "--single", extra={"FV_BIN_SHA256": "0" * 64})
    assert r.returncode != 0 and "not the pinned binary" in r.stderr
    assert not _calls(fv)


def test_delete_gnomad_src_after_success(fv):
    g = fv["src"] / "gnomad.exomes.v4.1.sites.chr21.vcf.bgz"
    r = _run(BUILD, fv, "--chrom", "21", "--delete-gnomad-src")
    assert r.returncode == 0, r.stderr
    assert not g.exists() and not Path(str(g) + ".tbi").exists()
    g22 = fv["src"] / "gnomad.exomes.v4.1.sites.chr22.vcf.bgz"
    r = _run(BUILD, fv, "--chrom", "22", "--delete-gnomad-src", extra={"FAIL_SOURCE": "phylop"})
    assert r.returncode != 0 and g22.exists()


def test_lsf_array_index_maps_to_chromosome(fv, tmp_path):
    env_sh = tmp_path / "env.sh"
    env_sh.write_text("".join(f'export {k}="{v}"\n' for k, v in fv["env"].items()
                              if k.startswith("FV_")))
    for idx, chrom in (("23", "X"), ("24", "Y"), ("7", "7")):
        r = subprocess.run(["bash", str(LSF)], capture_output=True, text=True,
                           env={**fv["env"], "FV_ENV": str(env_sh), "LSB_JOBINDEX": idx})
        assert r.returncode == 0, r.stderr
        assert f"-> chr{chrom}" in r.stderr
        assert (fv["db"] / f"gnomad_chr{chrom}.osa2.done").exists()


def _build_all(fv):
    assert _run(BUILD, fv, "--single").returncode == 0
    for c in CHROMS:
        r = _run(BUILD, fv, "--chrom", c)
        assert r.returncode == 0, r.stderr


def test_check_full_then_quick(fv):
    _build_all(fv)
    r = _run(CHECK, fv)
    assert r.returncode == 0, r.stderr
    assert "exactly 100 databases" in r.stderr
    sha = fv["db"] / "SA_DB_SHA256.txt"
    assert len(sha.read_text().splitlines()) == 100
    assert r.stdout.strip() == f"sa_db_sha256 {_sha(sha)}"
    man = (fv["db"] / "SA_DB_MANIFEST.tsv").read_text().splitlines()
    assert man[0].split("\t")[:3] == ["file", "size", "sha256"] and len(man) == 101
    q = _run(CHECK, fv, "--quick")
    assert q.returncode == 0, q.stderr
    assert q.stdout.strip() == r.stdout.strip()


def test_check_rejects_unexpected_and_missing(fv):
    _build_all(fv)
    assert _run(CHECK, fv).returncode == 0
    stray = fv["db"] / "gnomad_genomes_chr1.osa2"
    stray.write_text("x")
    r = _run(CHECK, fv, "--quick")
    assert r.returncode != 0 and "unexpected database file" in r.stderr
    stray.unlink()
    (fv["db"] / "revel_chr5.osa2.done").unlink()
    r = _run(CHECK, fv, "--quick")
    assert r.returncode != 0 and "revel_chr5.osa2" in r.stderr


def test_check_quick_detects_changed_db(fv):
    _build_all(fv)
    assert _run(CHECK, fv).returncode == 0
    (fv["db"] / "clinvar.osa2").write_text("a different database\n")
    r = _run(CHECK, fv, "--quick")
    assert r.returncode != 0 and "clinvar.osa2" in r.stderr


def test_check_rejects_non_exome_gnomad(fv):
    _build_all(fv)
    d = fv["db"] / "gnomad_chr3.osa2.done"
    d.write_text(d.read_text().replace("gnomad.exomes.v4.1", "gnomad.joint.v4.1"))
    r = _run(CHECK, fv, "--quick")
    assert r.returncode != 0 and "not built from gnomAD v4.1 exomes" in r.stderr


def test_check_fails_when_stack_or_omim_fails(fv):
    _build_all(fv)
    r = _run(CHECK, fv, extra={"STACK_FAIL": "1"})
    assert r.returncode != 0 and not (fv["db"] / "SA_DB_SHA256.txt").exists()
    r = _run(CHECK, fv, extra={"STACK_NO_OMIM": "1"})
    assert r.returncode != 0 and "omim" in r.stderr
    r = _run(CHECK, fv, "--quick")
    assert r.returncode != 0 and "run check_sa_db.sh (full) once" in r.stderr


@pytest.mark.parametrize("script", [BUILD, CHECK, LSF, HERE.parent / "fetch_sa_sources.sh",
                                    HERE.parent / "setup_fastvep.sh"])
def test_scripts_are_syntactically_valid(script):
    assert subprocess.run(["bash", "-n", str(script)]).returncode == 0


def _bump(path: Path, text: str):
    """Rewrite a source and move its mtime, as a new download would."""
    path.write_text(text)
    st = path.stat()
    os.utime(path, (st.st_atime, st.st_mtime + 100))


def test_changed_source_rebuilds_without_force(fv):
    assert _run(BUILD, fv, "--single").returncode == 0
    assert _run(BUILD, fv, "--chrom", "21").returncode == 0
    n = len(_calls(fv))
    assert _run(BUILD, fv, "--single").returncode == 0
    assert len(_calls(fv)) == n, "nothing changed, nothing should rebuild"

    _bump(fv["src"] / "clingen_gene_validity_20261110.csv", "a newer export\n")
    _bump(fv["src"] / "revel" / "revel_chr21.csv", "h\nanother row\n")
    r = _run(BUILD, fv, "--single")
    assert r.returncode == 0 and "[rebuild] omim.oga: source changed" in r.stderr
    r2 = _run(BUILD, fv, "--chrom", "21")
    assert r2.returncode == 0 and "[rebuild] revel_chr21.osa2: source changed" in r2.stderr
    new = [c.split("--source ")[1].split()[0] for c in _calls(fv)[n:]]
    assert sorted(new) == ["omim", "revel"], new


def test_use_omim_takes_effect_without_force(fv):
    (fv["src"] / "genemap2.txt").write_text("data\n")
    assert _run(BUILD, fv, "--single").returncode == 0
    n = len(_calls(fv))
    r = _run(BUILD, fv, "--single", "--use-omim")
    assert r.returncode == 0, r.stderr
    new = _calls(fv)[n:]
    assert len(new) == 1 and "--source omim" in new[0] and "genemap2.txt" in new[0]
    assert "[rebuild] omim.oga: source file was" in r.stderr
    # and back to ClinGen
    r = _run(BUILD, fv, "--single")
    assert r.returncode == 0 and "derived/clingen_gdv.tsv" in _calls(fv)[-1]


def test_deleted_gnomad_source_still_counts_as_current(fv):
    assert _run(BUILD, fv, "--chrom", "21", "--delete-gnomad-src").returncode == 0
    n = len(_calls(fv))
    r = _run(BUILD, fv, "--chrom", "21")
    assert r.returncode == 0, r.stderr
    assert len(_calls(fv)) == n


def test_new_binary_rebuilds(fv, tmp_path):
    assert _run(BUILD, fv, "--single").returncode == 0
    b2 = tmp_path / "fastvep2"
    b2.write_text(Path(fv["env"]["FV_BIN"]).read_text() + "# another build\n")
    b2.chmod(0o755)
    n = len(_calls(fv))
    r = _run(BUILD, fv, "--single", extra={"FV_BIN": str(b2), "FV_BIN_SHA256": _sha(b2)})
    assert r.returncode == 0 and "built by another fastvep" in r.stderr
    assert len(_calls(fv)) - n == 5
