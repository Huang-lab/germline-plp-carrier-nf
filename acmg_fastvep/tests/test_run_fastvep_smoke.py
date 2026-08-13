"""Smoke test: drive fastVEP end-to-end through run_fastvep.sh (plain mode).

Skipped automatically when `fastvep` is not on PATH, so `pytest` stays green on
machines without the binary (e.g. CI). When fastvep IS present, it annotates the
fastVEP repo's own shipped fixture (tests/test.vcf + tests/test.gff3, BRCA1/TP53
on chr17) and asserts a real CSQ header + CSQ-annotated records come out — i.e.
proves the wrapper actually runs the tool.

Point FASTVEP_SRC at a fastVEP checkout if it is not at /workspace/fastvep.
"""
from __future__ import annotations
import os
import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
RUNNER = HERE.parent / "run_fastvep.sh"
FASTVEP_SRC = Path(os.environ.get("FASTVEP_SRC", "/workspace/fastvep"))
TEST_VCF = FASTVEP_SRC / "tests" / "test.vcf"
TEST_GFF3 = FASTVEP_SRC / "tests" / "test.gff3"

_have_fastvep = shutil.which("fastvep") is not None
_have_fixture = TEST_VCF.is_file() and TEST_GFF3.is_file()


@pytest.mark.skipif(not _have_fastvep, reason="fastvep binary not on PATH")
@pytest.mark.skipif(not _have_fixture, reason=f"fastVEP fixture not found under {FASTVEP_SRC}")
def test_run_fastvep_plain(tmp_path):
    out = tmp_path / "out"
    r = subprocess.run(
        ["bash", str(RUNNER), "-i", str(TEST_VCF), "-o", str(out),
         "--gff3", str(TEST_GFF3)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, f"run_fastvep.sh failed:\nSTDOUT:{r.stdout}\nSTDERR:{r.stderr}"
    import gzip
    vcf = out / "test.fastvep.vcf.gz"
    assert vcf.is_file(), f"expected {vcf}; dir={list(out.iterdir()) if out.exists() else 'MISSING'}"
    text = gzip.open(vcf, "rt", encoding="utf-8").read()
    assert "##INFO=<ID=CSQ," in text, "no CSQ header in fastVEP output"
    # at least one non-comment record carries a CSQ= annotation
    assert any(("CSQ=" in ln) for ln in text.splitlines() if not ln.startswith("#")), \
        "no CSQ-annotated records in fastVEP output"


def test_runner_exists_and_parses_help():
    """Always-on: the wrapper exists and its usage/guards are wired (exit 2 on no args)."""
    assert RUNNER.is_file()
    r = subprocess.run(["bash", str(RUNNER)], capture_output=True, text=True)
    assert r.returncode == 2, "expected usage/exit-2 when required args are missing"
    assert "Usage: run_fastvep.sh" in r.stderr
