#!/usr/bin/env bash
# ClinVar-only path on a single-chromosome (chr21) synthetic fixture.
# Emulates: NORM_QC -> VEP -> VALIDATE_CHUNK (classifiers=clinvar)
#           -> CLINVAR_CLASSIFY -> PER_GENE_QC + CARRIER_MATRIX -> CONCAT.
#
# On Minerva:  nextflow run . -profile minerva -params-file params/msm.yaml \
#              --classifiers clinvar --input_vcfs '.../chr21.*.vcf.gz'
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export PYTHONPATH="$ROOT"
WORK="$(mktemp -d)"
cd "$WORK"

CHUNK="chr21_chunk"
IN="$ROOT/tests/synthetic/chunks_chr21/${CHUNK}.vcf"

# NORM_QC (fallback: no bcftools).
cp "$IN" "${CHUNK}.norm.vcf"

# VEP_ANNOTATE (fallback: synthetic_vep.py).
"$ROOT/bin/synthetic_vep.py" --in "${CHUNK}.norm.vcf" --out "${CHUNK}.vep.vcf"

# VALIDATE_CHUNK with classifiers=clinvar → only requires ClinVar CSQ subfields.
"$ROOT/bin/validate_chunk.py" \
    --vcf "${CHUNK}.vep.vcf" \
    --min-symbol-fraction 0.5 \
    --classifiers clinvar > "${CHUNK}.validate.json"
grep -q '"ok": true' "${CHUNK}.validate.json" || { cat "${CHUNK}.validate.json"; exit 1; }

# CLINVAR_CLASSIFY (min-stars=2).
"$ROOT/bin/classify_clinvar.py" \
    --vcf "${CHUNK}.vep.vcf" \
    --min-stars 2 \
    --out "${CHUNK}.clinvar_plp.tsv"

# CARRIER_MATRIX: build BED, extract GTs, matrix.
python3 - <<PY
import csv
pos=set()
for row in csv.DictReader(open("${CHUNK}.clinvar_plp.tsv"), delimiter="\t"):
    if int(row.get("is_clinvar_PLP","0") or 0):
        pos.add((row["chr"], row["pos"]))
with open("qualifying.pos.txt","w") as out:
    for c,p in sorted(pos):
        out.write(f"{c}\t{p}\n")
PY
"$ROOT/bin/fixture_gt_extract.py" --vcf "${CHUNK}.vep.vcf" --positions qualifying.pos.txt --out gt.tsv

# ACMG and AM inputs OMITTED (--acmg/--am unset).
"$ROOT/bin/build_carrier_matrix.py" \
    --gt gt.tsv \
    --clinvar "${CHUNK}.clinvar_plp.tsv" \
    --keep "$ROOT/tests/synthetic/samples.keep.chr21.txt" \
    --out-long "${CHUNK}.carrier_matrix.tsv"

# CONCAT (single chunk here, but exercise the code path).
head -n1 "${CHUNK}.carrier_matrix.tsv" > carrier_matrix.tsv
tail -n +2 "${CHUNK}.carrier_matrix.tsv" >> carrier_matrix.tsv

# MANIFEST.
"$ROOT/bin/make_manifest.py" \
    --pipeline-sha SMOKE \
    --classifiers clinvar \
    --clinvar-release SYNTHETIC \
    --out manifest.json

echo "===== validate.json ====="
cat "${CHUNK}.validate.json"; echo
echo "===== clinvar_plp.tsv ====="
cat "${CHUNK}.clinvar_plp.tsv"
echo "===== carrier_matrix.tsv ====="
cat carrier_matrix.tsv
echo "===== manifest.json ====="
cat manifest.json; echo
echo "CLINVAR-ONLY SMOKE OK"
