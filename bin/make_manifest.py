#!/usr/bin/env python3
"""Emit results/manifest.json capturing versions + provenance for this run."""
from __future__ import annotations
import argparse
import json
import os
import sys


def _safe_read(path: str) -> str:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return fh.read().strip()
    except OSError:
        return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--versions-file", default="", help="resources/versions.txt with resource release info")
    ap.add_argument("--pipeline-sha", default="")
    ap.add_argument("--reference-build", default="GRCh38")
    ap.add_argument("--vep-cache", default="")
    ap.add_argument("--vep-version", default="")
    ap.add_argument("--clinvar-release", default="")
    ap.add_argument("--gnomad-version", default="")
    ap.add_argument("--am-data-version", default="")
    ap.add_argument("--am-calibration-version", default="")
    ap.add_argument("--container-digests", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    manifest = {
        "reference_build": args.reference_build,
        "vep": {"version": args.vep_version, "cache": args.vep_cache},
        "clinvar_release": args.clinvar_release,
        "gnomad_version": args.gnomad_version,
        "alphamissense": {
            "data_version": args.am_data_version,
            "calibration_version": args.am_calibration_version,
        },
        "container_digests": args.container_digests,
        "pipeline_git_sha": args.pipeline_sha,
        "resource_versions_raw": _safe_read(args.versions_file),
    }
    with open(args.out, "w", encoding="utf-8") as out:
        json.dump(manifest, out, indent=2)
    print(args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
