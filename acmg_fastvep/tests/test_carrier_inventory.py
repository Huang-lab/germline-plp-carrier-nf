"""carrier_inventory.py and the inventory check in run_acmg_carriers_chr.sh.

The driver used to build a carrier table from whatever *.fastvep.vcf.gz it found, exit 0 and write
DONE over a quietly smaller (or mixed-binary) cohort. These tests build chromosome directories the
way run_fastvep_batch.sh leaves them (one test uses the real batch driver with a stub fastvep, so
the provenance and manifest formats cannot drift unnoticed) and check what the driver does when
the directory is not what the batch driver wrote.
"""
from __future__ import annotations

import importlib.util
import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ACMG = HERE.parent
DRIVER = ACMG / "run_acmg_carriers_chr.sh"
BATCH = ACMG / "run_fastvep_batch.sh"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py" if (HERE / f"{name}.py").exists() else ACMG / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pick = _load("test_fastvep_pick")          # write_vcf / GOOD: fastVEP v0.5.0-shaped CSQ fixtures
inv = _load("carrier_inventory")

PROV = """mode             acmg
fastvep_sha256   {sha}
sa_db_sha256     {bb}
input_sha256     {inp}
transcript_cache_sha256 {cc}
fasta            /ref/GRCh38.fa
acmg_config_sha256 none
fastvep_flags    --output-format vcf --acmg --flag-pick-allele-gene
chunk            {n}
"""


def make_dir(tmp_path, names=("c1", "c2", "c3"), *, manifest=True, run_prov=True, prov=True,
             status=None, bins=None, in_dir=True):
    d = tmp_path / "results-fastvep-chr21_fv050"
    d.mkdir()
    inp = tmp_path / "norm_qc"
    inp.mkdir()
    for n in names:
        pick.write_vcf(d / f"{n}.fastvep.vcf.gz", pick.GOOD, gz=True)
        (inp / f"{n}.norm.vcf.gz").write_bytes(b"x")
        if prov:
            (d / f"{n}.fastvep.provenance").write_text(PROV.format(
                sha=(bins or {}).get(n, "aa" * 32), bb="bb" * 32, cc="cc" * 32, inp=n * 10, n=n))
    if manifest:
        (d / "batch_manifest.tsv").write_text("".join(f"{n}\t{(status or {}).get(n, 'ok')}\n" for n in names))
    if run_prov:
        (d / "fastvep_run_provenance.txt").write_text(f"mode                   acmg\nin_dir                 {inp}\n")
    if not in_dir:
        shutil.rmtree(inp)
    return d, inp


def drive(tmp_path, d, *extra):
    out = tmp_path / "out"
    res = subprocess.run(["bash", str(DRIVER), "--in-dir", str(d), "-o", str(out), "-j", "1", *extra],
                         capture_output=True, text=True)
    return res, out


def done_detail(out):
    return (out / "DONE").read_text() if (out / "DONE").exists() else ""


# ---------------------------------------------------------------- the module
def test_complete_directory_has_no_problem(tmp_path):
    d, _ = make_dir(tmp_path)
    assert inv.check(str(d)) == ([], [], "checked")


def test_manifest_chunk_without_a_vcf(tmp_path):
    d, _ = make_dir(tmp_path)
    (d / "c2.fastvep.vcf.gz").unlink()
    problems, _, _ = inv.check(str(d))
    assert any("batch_manifest.tsv lists 1 chunk(s) with no .fastvep.vcf.gz here: c2" in p for p in problems)


def test_missing_chunk_is_found_through_the_input_list_when_there_is_no_manifest(tmp_path):
    d, _ = make_dir(tmp_path, manifest=False)
    (d / "c2.fastvep.vcf.gz").unlink()
    problems, _, _ = inv.check(str(d))
    assert any("input VCF(s)" in p and "c2" in p for p in problems)


def test_failed_and_submitted_chunks(tmp_path):
    d, _ = make_dir(tmp_path, status={"c2": "FAILED", "c3": "submitted"})
    problems, _, _ = inv.check(str(d))
    assert any("c2 as 'FAILED'" in p for p in problems)
    assert any("c3 as 'submitted'" in p and "without --lsf" in p for p in problems)


def test_chunks_from_two_binaries(tmp_path):
    d, _ = make_dir(tmp_path, bins={"c2": "dd" * 32})
    problems, _, _ = inv.check(str(d))
    assert any("not all made with the same fastvep_sha256" in p for p in problems)
    assert not any("input_sha256" in p for p in problems)      # per chunk by design


def test_vcf_without_provenance(tmp_path):
    d, _ = make_dir(tmp_path)
    (d / "c2.fastvep.provenance").unlink()
    problems, _, _ = inv.check(str(d))
    assert any("1 chunk(s) have no .fastvep.provenance" in p and "c2" in p for p in problems)


def test_vcf_the_manifest_does_not_list(tmp_path):
    d, inp = make_dir(tmp_path)
    pick.write_vcf(d / "c9.fastvep.vcf.gz", pick.GOOD, gz=True)
    (d / "c9.fastvep.provenance").write_text((d / "c1.fastvep.provenance").read_text())
    problems, _, _ = inv.check(str(d))
    assert any("not in batch_manifest.tsv" in p and "c9" in p for p in problems)


def test_unreadable_in_dir_is_a_note_not_a_failure(tmp_path):
    d, _ = make_dir(tmp_path, in_dir=False)
    problems, notes, level = inv.check(str(d))
    assert problems == [] and level == "checked"
    assert any("not readable from here" in n for n in notes)


def test_directory_without_any_record_is_legacy_not_a_pass(tmp_path):
    d, _ = make_dir(tmp_path, manifest=False, run_prov=False, prov=False)
    problems, notes, level = inv.check(str(d))
    assert problems == [] and level == "legacy-unchecked"
    assert any("completeness cannot be checked" in n for n in notes)


# ---------------------------------------------------------------- the driver
def test_driver_complete_directory_records_chunk_count_in_done(tmp_path):
    d, _ = make_dir(tmp_path)
    res, out = drive(tmp_path, d)
    assert res.returncode == 0, res.stderr
    assert "chunks=3 inventory=checked" in done_detail(out)


def test_driver_records_mode_and_pipeline_commit_in_done(tmp_path):
    """Issue #9: the table says which code and mode made it."""
    d, _ = make_dir(tmp_path)
    res, out = drive(tmp_path, d)
    assert res.returncode == 0, res.stderr
    detail = done_detail(out)
    rev = subprocess.run(["git", "-C", str(DRIVER.parent), "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    assert "mode=fastvep_pick" in detail
    assert f"pipeline_rev={rev}" in detail and len(rev) == 40
    assert "dirty=yes" in detail or "dirty=no" in detail


def test_driver_refuses_a_missing_chunk_before_building_anything(tmp_path):
    d, _ = make_dir(tmp_path)
    (d / "c2.fastvep.vcf.gz").unlink()
    res, out = drive(tmp_path, d)
    assert res.returncode == 4, res.stderr
    assert "no .fastvep.vcf.gz here: c2" in res.stderr
    assert (out / "FAILED").exists() and not (out / "DONE").exists()
    assert "chunk inventory check failed" in (out / "FAILED").read_text()
    assert not (out / "acmg_carriers.tsv").exists()


def test_driver_refuses_a_failed_chunk_whose_old_output_is_still_there(tmp_path):
    d, _ = make_dir(tmp_path, status={"c2": "FAILED"})
    res, out = drive(tmp_path, d)
    assert res.returncode == 4 and "'FAILED'" in res.stderr
    assert (out / "FAILED").exists()


def test_driver_refuses_mixed_binaries(tmp_path):
    d, _ = make_dir(tmp_path, bins={"c1": "dd" * 32})
    res, out = drive(tmp_path, d)
    assert res.returncode == 4 and "same fastvep_sha256" in res.stderr


def test_driver_legacy_directory_still_builds_and_says_it_was_not_checked(tmp_path):
    d, _ = make_dir(tmp_path, manifest=False, run_prov=False, prov=False)
    res, out = drive(tmp_path, d)
    assert res.returncode == 0, res.stderr
    assert "completeness cannot be checked" in res.stderr
    assert "inventory=legacy-unchecked" in done_detail(out)


def test_driver_allow_unchecked_inputs_builds_and_records_the_override(tmp_path):
    d, _ = make_dir(tmp_path)
    (d / "c2.fastvep.vcf.gz").unlink()
    res, out = drive(tmp_path, d, "--allow-unchecked-inputs")
    assert res.returncode == 0, res.stderr
    assert "WARN: --allow-unchecked-inputs" in res.stderr
    assert "chunks=2 inventory=overridden" in done_detail(out)


@pytest.mark.parametrize("field", inv.SAME_ACROSS_CHUNKS)
def test_every_origin_field_is_compared_across_chunks(tmp_path, field):
    d, _ = make_dir(tmp_path)
    p = d / "c2.fastvep.provenance"
    p.write_text("".join(f"{k:<16} CHANGED\n" if k == field else f"{k:<16} {v}\n"
                         for k, v in inv.read_kv(str(p)).items()))
    problems, _, _ = inv.check(str(d))
    assert any(f"same {field}" in x for x in problems), problems


# ------------------------------------- the real batch driver writes the directory, a stub is fastvep
def _stub_fastvep(tmp_path, vcf):
    p = tmp_path / "fastvep"
    p.write_text('#!/usr/bin/env bash\n[ "$1" = --version ] && { echo "fastvep 0.5.0"; exit 0; }\ncat "%s"\n' % vcf)
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return p


def _real_batch(tmp_path):
    """Run run_fastvep_batch.sh (stub fastvep) over three chunks; returns (dir, rerun())."""
    fixture = pick.write_vcf(tmp_path / "fixture.vcf", pick.GOOD)
    fv = _stub_fastvep(tmp_path, fixture)
    for n in ("gff3", "fa", "cache"):
        (tmp_path / n).write_text("x\n")
    inp = tmp_path / "norm_qc"
    inp.mkdir()
    for n in ("c1", "c2", "c3"):
        (inp / f"{n}.norm.vcf.gz").write_bytes(b"x")
    d = tmp_path / "results-fastvep-chr21_fv050"
    sha = subprocess.run(["bash", "-c", f"(sha256sum {fv} 2>/dev/null || shasum -a 256 {fv}) | cut -d' ' -f1"],
                         capture_output=True, text=True).stdout.strip()
    env = {**os.environ, "FV_BIN": str(fv), "FV_BIN_SHA256": sha, "FV_GFF3": str(tmp_path / "gff3"),
           "FV_FASTA": str(tmp_path / "fa"), "FV_CACHE": str(tmp_path / "cache")}

    def rerun():
        b = subprocess.run(["bash", str(BATCH), "--in-dir", str(inp), "-o", str(d), "--hgvs"],
                           capture_output=True, text=True, env=env)
        assert b.returncode == 0, b.stderr

    rerun()
    return d, rerun


def test_directory_written_by_the_real_batch_driver_passes_and_then_fails_without_a_chunk(tmp_path):
    d, _ = _real_batch(tmp_path)
    assert sorted(p.name for p in d.glob("*.fastvep.provenance")) == [f"c{i}.fastvep.provenance" for i in (1, 2, 3)]
    res, out = drive(tmp_path, d)
    assert res.returncode == 0, res.stderr
    assert "chunks=3 inventory=checked" in done_detail(out)
    (d / "c2.fastvep.vcf.gz").unlink()
    res, out = drive(tmp_path, d)
    assert res.returncode == 4 and "c2" in res.stderr and (out / "FAILED").exists()


def test_directory_after_a_rerun_that_skipped_every_chunk_still_passes(tmp_path):
    d, rerun = _real_batch(tmp_path)
    rerun()                                         # all three chunks skipped: manifest says skipped-already-done
    assert (d / "batch_manifest.tsv").read_text().count("skipped-already-done") == 3
    res, out = drive(tmp_path, d)
    assert res.returncode == 0, res.stderr
    assert "chunks=3 inventory=checked" in done_detail(out)
