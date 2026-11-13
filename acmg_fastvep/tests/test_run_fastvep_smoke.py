"""Smoke test: drive fastVEP end-to-end through run_fastvep.sh (plain mode).

Skipped unless FV_ENV points at an env.sh written by setup_fastvep.sh (the
pinned binary, transcript cache and FASTA), so `pytest` stays green on
machines without them (e.g. CI, a laptop). With FV_ENV set (on Minerva), it
annotates the fastVEP repo's own shipped fixture (tests/test.vcf, BRCA1/TP53 on
chr17) and asserts a real CSQ header + CSQ-annotated records come out — i.e.
proves the wrapper actually runs the tool.

FASTVEP_SRC defaults to $FV_SRC from env.sh (the v0.5.0 checkout).
"""
from __future__ import annotations
import os
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
RUNNER = HERE.parent / "run_fastvep.sh"
FV_ENV = os.environ.get("FV_ENV", "")
FASTVEP_SRC = Path(os.environ.get("FASTVEP_SRC") or os.environ.get("FV_SRC") or "/workspace/fastvep")
TEST_VCF = FASTVEP_SRC / "tests" / "test.vcf"

_have_env = bool(FV_ENV) and Path(FV_ENV).is_file()
_have_fixture = TEST_VCF.is_file()


@pytest.mark.skipif(not _have_env, reason="FV_ENV (setup_fastvep.sh env.sh) not set")
@pytest.mark.skipif(not _have_fixture, reason=f"fastVEP fixture not found under {FASTVEP_SRC}")
def test_run_fastvep_plain(tmp_path):
    out = tmp_path / "out"
    r = subprocess.run(
        ["bash", str(RUNNER), "-i", str(TEST_VCF), "-o", str(out)],
        capture_output=True, text=True, env={**os.environ, "FV_ENV": FV_ENV},
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
