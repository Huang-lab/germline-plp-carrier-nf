"""Functional tests for build_sa_dbs.sh (supplementary-annotation DB build).

The contract under test is *reporting*, not building: fastVEP does not fail when
a configured SA source is absent, it just returns no annotation for it, so an
incomplete sa_db changes ACMG calls with nothing on stderr. This script is the
only place a human sees which of the nine sources actually exist, so every one
of them must report a status on every run.

fastvep is stubbed — the point is which sources the script accounts for, not
what sa-build produces.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "build_sa_dbs.sh"

# The nine sources the ACMG classifier reads.
SOURCES = [
    "clinvar", "revel", "gnomad_genes", "clinvar_protein",
    "omim", "repeatmasker", "phylop", "spliceai", "gnomad",
]

# filename in sa_src -> the source it feeds
SRC_FILES = {
    "revel_with_transcript_ids": "revel",
    "gnomad_v4_constraint_metrics.tsv": "gnomad_genes",
    "variant_summary.txt.gz": "clinvar_protein",
    "genemap2.txt": "omim",
    "repeatmasker.bed": "repeatmasker",
    "hg38.phyloP100way.wigFix.gz": "phylop",
    "spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz": "spliceai",
    "gnomad.exomes.v4.1.sites.chr21.vcf.bgz": "gnomad",
}

STUB = """#!/usr/bin/env bash
out=""
while [ $# -gt 0 ]; do case "$1" in -o) out="$2"; shift 2;; *) shift;; esac; done
[ -n "$out" ] && : > "$out.oga"
exit 0
"""


def _env(tmp_path: Path, *, populate: bool, clinvar: bool = False):
    root = tmp_path / "root"
    src = root / "sa_src"
    src.mkdir(parents=True)
    binp = tmp_path / "fastvep"
    binp.write_text(STUB)
    binp.chmod(0o755)

    env = {
        "PATH": "/usr/bin:/bin",
        "FASTVEP": str(binp),
        "FV_ROOT": str(root),
    }
    if populate:
        for name in SRC_FILES:
            (src / name).write_text("data\n")
    if clinvar:
        cv = tmp_path / "clinvar_20260101.vcf.gz"
        cv.write_text("data\n")
        env["CLINVAR"] = str(cv)
    return env


def _run(env, extra=None):
    return subprocess.run(
        ["bash", str(SCRIPT)],
        env={**env, **(extra or {})},
        capture_output=True, text=True,
    )


def test_every_source_reports_a_status_when_none_are_present(tmp_path):
    """The regression: an absent source used to short-circuit silently.

    On a fresh machine the script named exactly one of the nine (revel, the only
    unconditional build_once call) and exited 0 having built nothing.
    """
    res = _run(_env(tmp_path, populate=False))
    out = res.stdout + res.stderr
    unreported = [s for s in SOURCES if s not in out]
    assert not unreported, f"sources that never reported a status: {unreported}\n{out}"


def test_summary_names_the_missing_sources(tmp_path):
    res = _run(_env(tmp_path, populate=False))
    out = res.stdout + res.stderr
    summary = [ln for ln in out.splitlines() if ln.startswith("[summary] sources MISSING:")]
    assert summary, f"no MISSING summary line\n{out}"
    for s in SOURCES:
        assert s in summary[0], f"{s} absent from {summary[0]!r}"


def test_summary_warns_that_acmg_calls_will_differ(tmp_path):
    res = _run(_env(tmp_path, populate=False))
    assert "ACMG calls made against this sa_db will differ" in res.stdout + res.stderr


def test_empty_source_tree_exits_zero(tmp_path):
    """`set -u` + an empty `chrs` array aborted the run on bash <= 4.3."""
    res = _run(_env(tmp_path, populate=False))
    assert res.returncode == 0, f"rc={res.returncode}\n{res.stdout}\n{res.stderr}"
    assert "unbound variable" not in res.stderr


def test_strict_mode_exits_nonzero_when_sources_are_missing(tmp_path):
    res = _run(_env(tmp_path, populate=False), {"SA_STRICT": "1"})
    assert res.returncode == 4, f"rc={res.returncode}\n{res.stderr}"


def test_full_source_tree_builds_and_reports_present(tmp_path):
    res = _run(_env(tmp_path, populate=True, clinvar=True))
    out = res.stdout + res.stderr
    present = [ln for ln in out.splitlines() if ln.startswith("[summary] sources present:")]
    assert present, f"no present summary\n{out}"
    # spliceai is deliberately not built by this script (fastVEP#101 chunk width),
    # so it is the one source that stays missing even with a complete sa_src.
    for s in ["clinvar", "revel", "gnomad_genes", "clinvar_protein",
              "omim", "repeatmasker", "phylop"]:
        assert s in present[0], f"{s} not reported present: {present[0]!r}"
    assert "gnomad_chr21" in present[0]
    assert "spliceai" not in present[0]


def test_spliceai_reports_differently_when_its_source_exists(tmp_path):
    """Source present but unbuildable here is a different state from absent."""
    absent = _run(_env(tmp_path / "a", populate=False))
    present = _run(_env(tmp_path / "b", populate=True))
    assert "no source file found" in absent.stderr
    assert "needs the sa-build(--format osa)" in present.stderr
