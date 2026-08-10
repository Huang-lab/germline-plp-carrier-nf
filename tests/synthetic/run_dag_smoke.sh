#!/usr/bin/env bash
# Emulates the Nextflow DAG on the synthetic fixture, invoking exactly the
# scripts each process invokes. Runs in the authoring environment where the
# `nextflow` binary is not installable (proxy policy blocks the installer).
#
# On Minerva, run the real pipeline with:
#     nextflow run . -profile test
#
# This shell driver is a smoke test for the wiring — it should stay in sync
# with the module scripts by construction.
set -euo pipefail
export PYTHONPATH="$(cd "$(dirname "$0")/../.." && pwd)"
ROOT="$PYTHONPATH"
cd "$(mktemp -d)"

CHUNK="chunk_001"
IN="$ROOT/tests/synthetic/chunks/${CHUNK}.vcf"

# NORM_QC (fallback path: no bcftools available)
cp "$IN" "${CHUNK}.norm.vcf"

# VEP_ANNOTATE (fallback: synthetic_vep.py)
"$ROOT/bin/synthetic_vep.py" --in "${CHUNK}.norm.vcf" --out "${CHUNK}.vep.vcf"

# VALIDATE_CHUNK
"$ROOT/bin/validate_chunk.py" --vcf "${CHUNK}.vep.vcf" --min-symbol-fraction 0.5 > "${CHUNK}.validate.json"
grep -q '"ok": true' "${CHUNK}.validate.json"

# CLINVAR_CLASSIFY
"$ROOT/bin/classify_clinvar.py" --vcf "${CHUNK}.vep.vcf" --min-stars 2 --out "${CHUNK}.clinvar_plp.tsv"

# ALPHAMISSENSE_CLASSIFY
"$ROOT/bin/classify_alphamissense.py" \
    --vcf "${CHUNK}.vep.vcf" \
    --calibration "$ROOT/tests/synthetic/am_calibration.tsv" \
    --min-strength PP3_Moderate \
    --out "${CHUNK}.am_plp.tsv"

# ACMG_FASTVEP (fallback: synthetic ACMG stub — no fastvep / SA dbs in the test env)
"$ROOT/bin/synthetic_acmg_fastvep.py" --in "${CHUNK}.vep.vcf" --out "${CHUNK}.acmg_plp.tsv"
test -s "${CHUNK}.acmg_plp.tsv" && echo "ACMG (fastVEP stub) table present"

# PER_GENE_QC (inlined python from modules/local/per_gene_qc.nf)
python3 - <<PY
import csv
from collections import Counter
counts = Counter(); plps = Counter()
for path, flag in [("${CHUNK}.clinvar_plp.tsv", "is_clinvar_PLP"),
                   ("${CHUNK}.acmg_plp.tsv", "is_acmg_PLP"),
                   ("${CHUNK}.am_plp.tsv", "is_AM_PLP")]:
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = row.get("gene", "") or "."
            counts[g] += 1
            if int(row.get(flag, "0") or 0):
                plps[(g, flag)] += 1
with open("${CHUNK}.qc_per_gene.tsv", "w") as out:
    out.write("gene\tn_annotated\tn_clinvar_PLP\tn_acmg_PLP\tn_AM_PLP\n")
    for g in sorted(counts):
        out.write(f"{g}\t{counts[g]}\t{plps[(g,'is_clinvar_PLP')]}\t{plps[(g,'is_acmg_PLP')]}\t{plps[(g,'is_AM_PLP')]}\n")
PY

# CARRIER_MATRIX
python3 - <<PY
import csv
pos=set()
for p,f in [("${CHUNK}.clinvar_plp.tsv","is_clinvar_PLP"),
            ("${CHUNK}.acmg_plp.tsv","is_acmg_PLP"),
            ("${CHUNK}.am_plp.tsv","is_AM_PLP")]:
    for row in csv.DictReader(open(p), delimiter="\t"):
        if int(row.get(f,"0") or 0):
            pos.add((row["chr"], row["pos"]))
with open("qualifying.pos.txt","w") as out:
    for c,p in sorted(pos):
        out.write(f"{c}\t{p}\n")
PY
"$ROOT/bin/fixture_gt_extract.py" --vcf "${CHUNK}.vep.vcf" --positions qualifying.pos.txt --out gt.tsv
"$ROOT/bin/build_carrier_matrix.py" \
    --gt gt.tsv \
    --clinvar "${CHUNK}.clinvar_plp.tsv" \
    --acmg "${CHUNK}.acmg_plp.tsv" \
    --am "${CHUNK}.am_plp.tsv" \
    --keep "$ROOT/tests/synthetic/samples.keep.txt" \
    --out-long "${CHUNK}.carrier_matrix.tsv"

# MANIFEST
"$ROOT/bin/make_manifest.py" \
    --pipeline-sha "SMOKE" \
    --reference-build GRCh38 \
    --clinvar-release SYNTHETIC \
    --gnomad-version SYNTHETIC \
    --out manifest.json

echo "----- carrier_matrix.tsv -----"
cat "${CHUNK}.carrier_matrix.tsv"
echo "----- qc_per_gene.tsv -----"
cat "${CHUNK}.qc_per_gene.tsv"
echo "----- manifest.json -----"
cat manifest.json
echo "SMOKE OK"
