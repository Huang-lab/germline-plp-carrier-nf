"""submit_fv050_pipeline.sh --dry-run: the three steps per chromosome and their LSF chaining."""
from __future__ import annotations

import gzip
import hashlib
import os
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ACMG = HERE.parent


def _world(tmp_path: Path):
    # A copy of the scripts whose check_sa_db.sh always passes (the real one needs
    # the 100-database set), so the batch's --acmg pre-flight can be exercised.
    scripts = tmp_path / "repo" / "acmg_fastvep"
    shutil.copytree(ACMG, scripts, ignore=shutil.ignore_patterns("tests", "__pycache__"))
    shutil.copytree(ACMG.parent / "plp_rules", tmp_path / "repo" / "plp_rules",
                    ignore=shutil.ignore_patterns("__pycache__"))
    (scripts / "check_sa_db.sh").write_text("#!/usr/bin/env bash\nexit 0\n")
    (scripts / "check_sa_db.sh").chmod(0o755)

    fake = tmp_path / "fastvep"
    fake.write_text('#!/usr/bin/env bash\n[ "$1" = --version ] && { echo "fastvep 0.5.0"; exit 0; }\nexit 0\n')
    fake.chmod(0o755)
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    (bin_ / "tabix").write_text("#!/usr/bin/env bash\nexit 0\n")
    (bin_ / "tabix").chmod(0o755)

    B = tmp_path / "batch"
    for c, chunks in (("1", 2), ("13", 1), ("X", 1)):
        d = B / f"results-chr{c}" / "norm_qc"
        d.mkdir(parents=True)
        for i in range(chunks):
            with gzip.open(d / f"s.{i}_chr{c}_1_100.norm.vcf.gz", "wt") as fh:
                fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    sa = tmp_path / "sa_db"
    sa.mkdir()
    (sa / "SA_DB_SHA256.txt").write_text("abc  clinvar.osa2\n")
    for name in ("g.gff3", "g.fa", "g.cache"):
        (tmp_path / name).write_text("x\n")
    envf = tmp_path / "env.sh"
    envf.write_text("".join(f"export {k}={v}\n" for k, v in {
        "FV_BIN": fake, "FV_BIN_SHA256": hashlib.sha256(fake.read_bytes()).hexdigest(),
        "FV_GFF3": tmp_path / "g.gff3", "FV_FASTA": tmp_path / "g.fa",
        "FV_CACHE": tmp_path / "g.cache", "FV_SA": sa}.items()))
    env = {k: v for k, v in os.environ.items() if not k.startswith("FV_")}
    env.update({"FV_ENV": str(envf), "MINERVA_ALLOCATION": "acc_test",
                "PATH": f"{bin_}:{env['PATH']}"})
    return scripts, B, env


def _submit(scripts, B, env, tmp_path, *extra):
    return subprocess.run(
        ["bash", str(scripts / "submit_fv050_pipeline.sh"), "--batch-dir", str(B),
         "--chroms", "1 13 17 X", "--skip", "13", "--run-tag", "T", "--dry-run",
         "--fv-out-root", str(tmp_path / "scratch"), "--log-dir", str(tmp_path / "logs"), *extra],
        env=env, capture_output=True, text=True, cwd=tmp_path)


def test_dry_run_chains_chunks_verify_carriers(tmp_path):
    scripts, B, env = _world(tmp_path)
    r = _submit(scripts, B, env, tmp_path)
    assert r.returncode == 0, r.stderr
    err = r.stderr
    bsubs = [l for l in (r.stdout + r.stderr).splitlines() if l.startswith("[dry-run] bsub")]
    # chr1: 2 chunk jobs + verify + carriers; chrX: 1 + 2; chr13 skipped; chr17 has no input
    assert len(bsubs) == 7, "\n".join(bsubs)
    assert sum("-J fv050-T-acmg-chr1-s." in l for l in bsubs) == 2
    assert sum("-J fv050-T-acmg-chrX-s." in l for l in bsubs) == 1
    assert not any("chr13" in l for l in bsubs)
    v1 = next(l for l in bsubs if "-J fv050-T-verify-chr1 " in l)
    assert "ended\\(fv050-T-acmg-chr1-\\*\\)" in v1 or "ended(fv050-T-acmg-chr1-*)" in v1
    assert "--lsf" not in v1 and "--acmg" in v1 and "module\\ load" in v1
    c1 = next(l for l in bsubs if "-J fv050-T-carriers-chr1 " in l)
    assert "done(fv050-T-verify-chr1)" in c1.replace("\\", "")
    assert "--mode\\ fastvep_pick" in c1 or "--mode fastvep_pick" in c1
    chunk = next(l for l in bsubs if "fv050-T-acmg-chr1-s." in l)
    # --mem 16000 over 4 threads: Minerva's rusage[mem] is per core
    assert "rusage[mem=4000]" in chunk.replace("\\", "")
    assert "rusage[mem=4000]" in v1.replace("\\", "") and "-n 8" in v1
    assert str(tmp_path / "scratch" / "results-fastvep-chr1_fv050") in chunk.replace("\\", "")
    assert "chr13: skipped (--skip)" in err and "no input:  17" in err


def test_nothing_left_means_no_wait_on_chunks(tmp_path):
    """A chromosome whose chunks are all done submits verify with no -w on chunk jobs."""
    scripts, B, env = _world(tmp_path)
    # make every chunk "done": run the batch once for real (no --lsf) with the fake binary
    for c in ("1", "X"):
        r = subprocess.run(["bash", str(scripts / "run_fastvep_batch.sh"), "--in-dir",
                            str(B / f"results-chr{c}" / "norm_qc"), "-o",
                            str(tmp_path / "scratch" / f"results-fastvep-chr{c}_fv050"),
                            "--hgvs", "--acmg"], env={**env, **_envfile(env)},
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
    r = _submit(scripts, B, env, tmp_path)
    assert r.returncode == 0, r.stderr
    bsubs = [l for l in (r.stdout + r.stderr).splitlines() if l.startswith("[dry-run] bsub")]
    assert not any("-acmg-chr" in l and "-J fv050-T-acmg" in l for l in bsubs)
    v1 = next(l for l in bsubs if "-J fv050-T-verify-chr1 " in l)
    assert "ended" not in v1


def _envfile(env):
    out = {}
    for line in Path(env["FV_ENV"]).read_text().splitlines():
        k, v = line.removeprefix("export ").split("=", 1)
        out[k] = v
    return out


def test_broken_python_refused_before_any_submission(tmp_path):
    scripts, B, env = _world(tmp_path)
    bad = tmp_path / "badbin"
    bad.mkdir()
    (bad / "python3").write_text("#!/usr/bin/env bash\nexit 127\n")
    (bad / "python3").chmod(0o755)
    r = _submit(scripts, B, {**env, "PATH": f"{bad}:{env['PATH']}"}, tmp_path)
    assert r.returncode == 3 and "python3" in r.stderr
    assert "[dry-run] bsub" not in r.stderr


def _with_bjobs(env, tmp_path, active: str, known: str):
    d = tmp_path / "lsfbin"
    d.mkdir(exist_ok=True)
    (d / "bjobs").write_text(
        "#!/usr/bin/env bash\n"
        f'case " $* " in *" -a "*) printf %s {known!r} ;; *) printf %s {active!r} ;; esac\n')
    (d / "bjobs").chmod(0o755)
    return {**env, "PATH": f"{d}:{env['PATH']}"}


def test_chromosome_with_queued_jobs_is_not_resubmitted(tmp_path):
    scripts, B, env = _world(tmp_path)
    env = _with_bjobs(env, tmp_path, "fv050-11161200-acmg-chr1-s.0_chr1_1_100\nfv050-acmg-chr10\n", "")
    r = _submit(scripts, B, env, tmp_path)
    assert r.returncode == 0, r.stderr
    bsubs = [l for l in (r.stdout + r.stderr).splitlines() if l.startswith("[dry-run] bsub")]
    assert not any("-chr1-" in l or "verify-chr1 " in l for l in bsubs)
    assert any("verify-chrX" in l for l in bsubs)
    assert "busy (earlier jobs still queued), not submitted: 1" in r.stderr


def test_reused_run_tag_is_refused(tmp_path):
    scripts, B, env = _world(tmp_path)
    env = _with_bjobs(env, tmp_path, "", "fv050-T-verify-chr1\n")
    r = _submit(scripts, B, env, tmp_path)
    assert r.returncode == 3 and "--run-tag" in r.stderr
