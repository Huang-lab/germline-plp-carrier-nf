"""DONE / FAILED markers in a carrier-table directory, and the check every reader makes.

A carrier run that finds a problem (PICK anomalies, a missing chunk) still
writes its acmg_variants.tsv / acmg_carriers.tsv, so the rows can be looked at,
and then exits non-zero. An exit code is easy to miss once the job has scrolled
by, while the tables look complete. So the directory itself says it:

  FAILED   written on a non-zero exit: reason, counts, exit code, UTC time
  DONE     written on success; a stale FAILED (or DONE) is removed first

Readers (extract_gene_carriers.py, compare_acmg_pick_vs_mane.py) refuse a dir
holding FAILED unless told --allow-pick-anomalies / --allow-failed, and warn
about a dir with neither marker (written before markers existed, or killed).
"""
from __future__ import annotations

import os
import sys
import time

DONE = "DONE"
FAILED = "FAILED"


def clear(d: str) -> None:
    for m in (DONE, FAILED):
        try:
            os.remove(os.path.join(d, m))
        except FileNotFoundError:
            pass


def _write(path: str, lines: list[str]) -> None:
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    os.replace(tmp, path)


def write_done(d: str, detail: str = "") -> None:
    clear(d)
    _write(os.path.join(d, DONE), [f"utc\t{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
                                   f"detail\t{detail}"])


def write_failed(d: str, reason: str, exit_code: int, counts: str = "") -> None:
    clear(d)
    _write(os.path.join(d, FAILED), [f"utc\t{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
                                     f"exit_code\t{exit_code}", f"reason\t{reason}",
                                     f"counts\t{counts}"])


def check_dir(d: str, allow_failed: bool, who: str) -> None:
    """Exit 4 if d holds FAILED (unless allow_failed); warn if it has no DONE."""
    failed = os.path.join(d, FAILED)
    if os.path.exists(failed):
        with open(failed, encoding="utf-8") as fh:
            body = fh.read().strip().replace("\n", "; ")
        if not allow_failed:
            print(f"ERROR: {who}: {d} holds a FAILED marker ({body}). Its tables were written by a "
                  "run that failed; refusing to read them. Fix and rerun, or pass "
                  "--allow-pick-anomalies to read them anyway.", file=sys.stderr)
            sys.exit(4)
        print(f"WARN: {who}: reading {d} despite FAILED ({body})", file=sys.stderr)
    elif not os.path.exists(os.path.join(d, DONE)):
        print(f"WARN: {who}: {d} has no DONE marker (written before markers existed, or the "
              "run was killed); its tables may be incomplete", file=sys.stderr)


if __name__ == "__main__":
    # Shell entry point for run_acmg_carriers_chr.sh:
    #   carrier_markers.py done <dir> [detail]
    #   carrier_markers.py failed <dir> <exit code> <reason> [counts]
    #   carrier_markers.py clear <dir>
    cmd, d = sys.argv[1], sys.argv[2]
    if cmd == "done":
        write_done(d, sys.argv[3] if len(sys.argv) > 3 else "")
    elif cmd == "failed":
        write_failed(d, sys.argv[4], int(sys.argv[3]), sys.argv[5] if len(sys.argv) > 5 else "")
    elif cmd == "clear":
        clear(d)
    else:
        sys.exit(f"unknown command {cmd}")
