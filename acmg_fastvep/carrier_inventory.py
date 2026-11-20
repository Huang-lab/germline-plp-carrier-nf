"""Is this chromosome directory of fastVEP chunks complete and of one origin?

run_acmg_carriers_chr.sh builds a carrier table from whatever *.fastvep.vcf.gz it finds in a
directory. A chunk whose job never ran, was purged, or failed on a rerun leaves the other
chunks behind, and the table is then well formed, exits 0 and gets a DONE marker over a
quietly smaller cohort. This module compares the directory with what the batch driver
(run_fastvep_batch.sh) recorded when it wrote it:

  <name>.fastvep.provenance   per chunk: binary, SA database, transcript cache, FASTA, flags, input hash
  batch_manifest.tsv          per chunk: ok | skipped-already-done | FAILED | submitted
  fastvep_run_provenance.txt  per run: in_dir, the directory of QC'd VCFs the chunks were made from

and reports every disagreement:
  * a chunk the manifest lists that has no VCF here, or a VCF the manifest does not list
  * a chunk the manifest says is FAILED or only `submitted` (its job's outcome was never recorded)
  * a VCF without a .fastvep.provenance (the batch driver removes it first on every rerun)
  * chunks made with different binaries / SA databases / caches / FASTAs / flags / modes
  * an input VCF in in_dir with no output here, or an output whose input is gone (when in_dir is readable)

A directory with none of those files (made before the batch driver recorded them) cannot be
checked; that is reported as `legacy-unchecked`, not as a pass.

CLI: carrier_inventory.py <chunk dir>   -> prints `checked` or `legacy-unchecked` on stdout;
     problems and notes on stderr; exit 4 on any problem.
"""
from __future__ import annotations

import glob
import os
import sys

OK_STATUS = ("ok", "skipped-already-done")
# every field of run_fastvep.sh's provenance_key() except input_sha256 (which is per chunk by design)
SAME_ACROSS_CHUNKS = ("mode", "fastvep_sha256", "sa_db_sha256", "transcript_cache_sha256", "fasta",
                      "acmg_config_sha256", "fastvep_flags")
VCF_SUFFIX = ".fastvep.vcf.gz"


def chunk_name(path: str) -> str:
    """The batch driver's chunk_name(): basename without .gz, .vcf, then .norm."""
    n = os.path.basename(path)
    for suf in (".gz", ".vcf", ".norm"):
        if n.endswith(suf):
            n = n[: -len(suf)]
    return n


def read_kv(path: str) -> dict:
    out: dict = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            parts = ln.rstrip("\n").split(None, 1)
            if parts and parts[0] not in out:
                out[parts[0]] = parts[1].strip() if len(parts) > 1 else ""
    return out


def read_manifest(path: str) -> dict:
    out: dict = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f and f[0]:
                out[f[0]] = f[1] if len(f) > 1 else ""
    return out


def input_names(in_dir: str) -> list:
    """The batch driver's input choice: *.norm.vcf.gz, else *.vcf.gz, else *.vcf."""
    for pat in ("*.norm.vcf.gz", "*.vcf.gz", "*.vcf"):
        hits = sorted(glob.glob(os.path.join(in_dir, pat)))
        if hits:
            return [chunk_name(h) for h in hits]
    return []


def _short(names, limit=8) -> str:
    names = sorted(names)
    return ", ".join(names[:limit]) + (f", ... ({len(names)} in all)" if len(names) > limit else "")


def check(d: str):
    """-> (problems, notes, level). level is 'checked' or 'legacy-unchecked'."""
    vcfs = sorted(glob.glob(os.path.join(d, "*" + VCF_SUFFIX)))
    names = {os.path.basename(v)[: -len(VCF_SUFFIX)] for v in vcfs}
    man_p = os.path.join(d, "batch_manifest.tsv")
    run_p = os.path.join(d, "fastvep_run_provenance.txt")
    prov = {n: os.path.join(d, n + ".fastvep.provenance") for n in names}
    have_prov = {n for n, p in prov.items() if os.path.exists(p)}
    if not os.path.exists(man_p) and not os.path.exists(run_p) and not have_prov:
        return [], ["no batch_manifest.tsv, fastvep_run_provenance.txt or .fastvep.provenance here "
                    "(made before the batch driver recorded them): completeness cannot be checked"], \
            "legacy-unchecked"
    problems: list = []
    notes: list = []

    if os.path.exists(man_p):
        man = read_manifest(man_p)
        bad = {n: s for n, s in man.items() if s not in OK_STATUS}
        for n, s in sorted(bad.items()):
            hint = ("; its job's outcome was never recorded: re-run run_fastvep_batch.sh without --lsf "
                    "(it skips finished chunks and rewrites the manifest)" if s == "submitted" else "")
            problems.append(f"batch_manifest.tsv lists chunk {n} as {s!r}{hint}")
        gone = sorted(set(man) - names)
        if gone:
            problems.append(f"batch_manifest.tsv lists {len(gone)} chunk(s) with no {VCF_SUFFIX} here: {_short(gone)}")
        extra = sorted(names - set(man))
        if extra:
            problems.append(f"{len(extra)} {VCF_SUFFIX} not in batch_manifest.tsv (output of a run whose "
                            f"manifest was replaced, or of a chunk no longer in the input): {_short(extra)}")
    else:
        notes.append("no batch_manifest.tsv: the chunk list cannot be compared with the batch driver's record")

    missing_prov = sorted(names - have_prov)
    if missing_prov and have_prov or (missing_prov and os.path.exists(man_p)):
        problems.append(f"{len(missing_prov)} chunk(s) have no .fastvep.provenance (the batch driver deletes it "
                        f"before a rerun and writes it last, so the output may be from an unfinished rerun): "
                        f"{_short(missing_prov)}")
    seen: dict = {}
    for n in sorted(have_prov):
        kv = read_kv(prov[n])
        for k in SAME_ACROSS_CHUNKS:
            seen.setdefault(k, {}).setdefault(kv.get(k, "<absent>"), []).append(n)
    for k, vals in seen.items():
        if len(vals) > 1:
            parts = "; ".join(f"{v[:16]}{'...' if len(v) > 16 else ''} in {len(ns)} chunk(s) ({_short(ns, 3)})"
                              for v, ns in sorted(vals.items(), key=lambda kv: -len(kv[1])))
            problems.append(f"chunks were not all made with the same {k}: {parts}")

    in_dir = ""
    if os.path.exists(run_p):
        in_dir = read_kv(run_p).get("in_dir", "")
    if in_dir and os.path.isdir(in_dir):
        want = set(input_names(in_dir))
        if not want:
            notes.append(f"in_dir {in_dir} holds no VCFs now: the input list cannot be compared")
        else:
            if want - names:
                problems.append(f"{len(want - names)} input VCF(s) in {in_dir} have no {VCF_SUFFIX} here: {_short(want - names)}")
            if names - want:
                problems.append(f"{len(names - want)} {VCF_SUFFIX} here have no input VCF in {in_dir}: {_short(names - want)}")
    elif in_dir:
        notes.append(f"in_dir {in_dir} (fastvep_run_provenance.txt) is not readable from here: the input list was not compared")
    else:
        notes.append("no in_dir recorded: the input list was not compared")
    return problems, notes, "checked"


def main(argv) -> int:
    if len(argv) != 2 or not os.path.isdir(argv[1]):
        print("usage: carrier_inventory.py <directory of *.fastvep.vcf.gz chunks>", file=sys.stderr)
        return 2
    problems, notes, level = check(argv[1])
    for n in notes:
        print(f"[inventory] note: {n}", file=sys.stderr)
    for p in problems:
        print(f"[inventory] PROBLEM: {p}", file=sys.stderr)
    if problems:
        print(f"ERROR: {argv[1]} is not a complete, single-origin set of chunks ({len(problems)} problem(s) above); "
              "refusing to build a carrier table that would silently under-count or mix the cohort. "
              "Fix the directory (re-run run_fastvep_batch.sh), or pass --allow-unchecked-inputs to build anyway.",
              file=sys.stderr)
        return 4
    print(level)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
