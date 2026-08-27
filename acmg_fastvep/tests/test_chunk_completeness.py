"""Regression tests for the two ways a chunk could go missing without a word.

Both failure modes were live before this file existed, and both produce a
well-formed artifact over a quietly smaller cohort — the shape of the batch3
dropped-chunk incident:

1. `run_fastvep.sh` streamed its output straight into the final path, so a chunk
   killed at the LSF walltime left a truncated `.fastvep.vcf.gz` behind. The
   truncation is a *valid* gzip stream (gzip closes cleanly when its stdin goes
   away), so `gzip -t` passes and only the record count is wrong.
   `run_fastvep_batch.sh` then read `[ -s ... ]` as "already done" and skipped
   the chunk on the resubmit that was meant to finish it.

2. `build_carrier_matrix.sh` globbed whatever classifier TSVs were on disk. The
   Nextflow path cannot do this — `carrier_gt.nf` joins each classifier TSV to
   its own chunk VCF on `chunk_id` — but the standalone script had no equivalent
   invariant, so a missing chunk just made the carrier matrix smaller.

Each test is written to fail against the pre-fix scripts. Only `bash` and
`python3` are needed; no fastvep, no bcftools.
"""
from __future__ import annotations

import gzip
import os
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ACMG = HERE.parent
RUNNER = ACMG / "run_fastvep.sh"
BATCH = ACMG / "run_fastvep_batch.sh"
CARRIER = ACMG / "build_carrier_matrix.sh"

ACMG_HEADER = (
    "chr\tpos\tref\talt\tgene\tacmg_label\tacmg_criteria\t"
    "n_pathogenic_criteria\tn_benign_criteria\tis_acmg_PLP"
)


def _fake_fastvep(tmp_path: Path, *, records: int, exit_code: int) -> Path:
    """A stand-in for the fastvep binary that emits `records` lines then exits.

    `exit_code != 0` reproduces a job killed partway through: the VCF it had
    written so far is a complete gzip member holding a prefix of the records.
    """
    p = tmp_path / "fake_fastvep"
    p.write_text(
        "#!/usr/bin/env bash\n"
        "# argv is the fastvep annotate command line; ignore it and emit a VCF.\n"
        'printf "##fileformat=VCFv4.2\\n"\n'
        'printf "##INFO=<ID=CSQ,Number=.,Type=String,Description=\\"x|ACMG|ACMG_CRITERIA\\">\\n"\n'
        'printf "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n"\n'
        f"for i in $(seq 1 {records}); do\n"
        '  printf "17\\t%s\\t.\\tA\\tG\\t.\\tPASS\\tCSQ=x|LP|PVS1\\n" "$i"\n'
        "done\n"
        f"exit {exit_code}\n"
    )
    p.chmod(0o755)
    return p


def _chunk_vcf(d: Path, name: str) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    vcf = d / f"{name}.norm.vcf.gz"
    with gzip.open(vcf, "wt", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n")
        fh.write("17\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n")
    return vcf


def test_killed_run_publishes_nothing(tmp_path):
    """A chunk that dies mid-write must leave no final output to be mistaken
    for a finished one."""
    indir = tmp_path / "norm_qc"
    outdir = tmp_path / "out"
    _chunk_vcf(indir, "chunk1")
    gff3 = tmp_path / "genes.gff3"
    gff3.write_text("##gff-version 3\n")
    fake = _fake_fastvep(tmp_path, records=200, exit_code=137)

    r = subprocess.run(
        ["bash", str(RUNNER), "-i", str(indir / "chunk1.norm.vcf.gz"),
         "-o", str(outdir), "--gff3", str(gff3), "--fastvep", str(fake)],
        capture_output=True, text=True,
    )
    assert r.returncode != 0, "a killed fastvep must fail the wrapper"

    published = outdir / "chunk1.fastvep.vcf.gz"
    assert not published.exists(), (
        f"{published.name} was published by a run that died: "
        f"{published.stat().st_size if published.exists() else 0} bytes. "
        "The batch driver reads its presence as 'already done' and will skip "
        "this chunk forever."
    )
    # And no partial left lying around to confuse the next run.
    assert not list(outdir.glob(".*partial")), "stale .partial left behind"


def test_truncated_output_is_a_valid_gzip_stream(tmp_path):
    """Documents why `[ -s ]` and `gzip -t` cannot be the completeness test.

    This is the property that made the bug invisible; it is not fixed, it is
    routed around by publishing on rename.
    """
    fake = _fake_fastvep(tmp_path, records=200, exit_code=137)
    out = tmp_path / "truncated.vcf.gz"
    with out.open("wb") as fh:
        producer = subprocess.Popen([str(fake)], stdout=subprocess.PIPE)
        gz = subprocess.Popen(["gzip", "-c"], stdin=producer.stdout, stdout=fh)
        producer.stdout.close()
        gz.wait()
        producer.wait()

    assert out.stat().st_size > 0
    assert subprocess.run(["gzip", "-t", str(out)]).returncode == 0, (
        "expected the truncated output to pass gzip -t"
    )


def test_batch_driver_reruns_a_chunk_left_partial(tmp_path):
    """After a killed chunk, re-running the batch driver must retry it."""
    indir = tmp_path / "norm_qc"
    outdir = tmp_path / "out"
    _chunk_vcf(indir, "chunk1")
    gff3 = tmp_path / "genes.gff3"
    gff3.write_text("##gff-version 3\n")
    env = dict(os.environ)

    # First attempt dies.
    dying = _fake_fastvep(tmp_path, records=50, exit_code=137)
    subprocess.run(
        ["bash", str(BATCH), "--in-dir", str(indir), "-o", str(outdir),
         "--gff3", str(gff3), "--fastvep", str(dying)],
        capture_output=True, text=True, env=env,
    )
    manifest = (outdir / "batch_manifest.tsv").read_text()
    assert "FAILED" in manifest, f"first attempt should record a failure, got:\n{manifest}"

    # Second attempt succeeds and must actually run, not skip.
    good = _fake_fastvep(tmp_path, records=50, exit_code=0)
    subprocess.run(
        ["bash", str(BATCH), "--in-dir", str(indir), "-o", str(outdir),
         "--gff3", str(gff3), "--fastvep", str(good)],
        capture_output=True, text=True, env=env,
    )
    manifest = (outdir / "batch_manifest.tsv").read_text()
    assert "skipped-already-done" not in manifest, (
        "the retry skipped a chunk whose only output came from a killed run:\n"
        + manifest
    )
    assert "\tok" in manifest, f"retry should have completed the chunk, got:\n{manifest}"


def test_carrier_matrix_refuses_a_missing_chunk(tmp_path):
    """One classifier TSV short of the chunk count must be an error, not a
    smaller matrix."""
    norm = tmp_path / "norm_qc"
    acmg = tmp_path / "acmg"
    acmg.mkdir()
    for name in ("chunk1", "chunk2", "chunk3"):
        _chunk_vcf(norm, name)
    # Only two of the three chunks produced a table.
    for name in ("chunk1", "chunk2"):
        (acmg / f"{name}.acmg_plp.tsv").write_text(
            ACMG_HEADER + "\n17\t100\tA\tG\tBRCA1\tPathogenic\tPVS1\t1\t0\t1\n"
        )

    r = subprocess.run(
        ["bash", str(CARRIER), "--norm-dir", str(norm),
         "--acmg-dir", str(acmg), "-o", str(tmp_path / "carrier_matrix.tsv")],
        capture_output=True, text=True,
    )
    assert r.returncode != 0, (
        "built a carrier matrix from 2 of 3 chunks and exited 0:\n"
        f"STDOUT:{r.stdout}\nSTDERR:{r.stderr}"
    )
    assert "3 chunk VCF" in r.stderr or "chunk table" in r.stderr, r.stderr
    assert not (tmp_path / "carrier_matrix.tsv").exists() or \
        (tmp_path / "carrier_matrix.tsv").stat().st_size == 0


def test_carrier_matrix_allows_a_missing_chunk_when_asked(tmp_path):
    """--allow-partial is the deliberate escape hatch, and it warns."""
    norm = tmp_path / "norm_qc"
    acmg = tmp_path / "acmg"
    acmg.mkdir()
    for name in ("chunk1", "chunk2"):
        _chunk_vcf(norm, name)
    (acmg / "chunk1.acmg_plp.tsv").write_text(
        ACMG_HEADER + "\n17\t100\tA\tG\tBRCA1\tPathogenic\tPVS1\t1\t0\t1\n"
    )

    r = subprocess.run(
        ["bash", str(CARRIER), "--norm-dir", str(norm), "--acmg-dir", str(acmg),
         "--allow-partial", "-o", str(tmp_path / "carrier_matrix.tsv")],
        capture_output=True, text=True,
    )
    assert "WARNING" in r.stderr, f"--allow-partial must still say so:\n{r.stderr}"


def test_carrier_matrix_refuses_mismatched_headers(tmp_path):
    """`tail -n +2` splices rows under the first chunk's header; disagreeing
    columns must stop the run rather than misalign silently."""
    norm = tmp_path / "norm_qc"
    acmg = tmp_path / "acmg"
    acmg.mkdir()
    for name in ("chunk1", "chunk2"):
        _chunk_vcf(norm, name)
    (acmg / "chunk1.acmg_plp.tsv").write_text(
        ACMG_HEADER + "\n17\t100\tA\tG\tBRCA1\tPathogenic\tPVS1\t1\t0\t1\n"
    )
    # Same columns, different order — the case a diff of the data would not show.
    reordered = "\t".join(
        ["chr", "pos", "ref", "alt", "gene", "acmg_label", "acmg_criteria",
         "n_benign_criteria", "n_pathogenic_criteria", "is_acmg_PLP"]
    )
    (acmg / "chunk2.acmg_plp.tsv").write_text(
        reordered + "\n17\t200\tC\tT\tBRCA2\tPathogenic\tPVS1\t0\t1\t1\n"
    )

    r = subprocess.run(
        ["bash", str(CARRIER), "--norm-dir", str(norm), "--acmg-dir", str(acmg),
         "-o", str(tmp_path / "carrier_matrix.tsv")],
        capture_output=True, text=True,
    )
    assert r.returncode != 0, "concatenated chunks whose columns disagree"
    assert "header mismatch" in r.stderr, r.stderr


@pytest.mark.parametrize("script", [RUNNER, BATCH, CARRIER])
def test_scripts_are_syntactically_valid(script):
    assert subprocess.run(["bash", "-n", str(script)]).returncode == 0
