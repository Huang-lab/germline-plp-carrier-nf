#!/usr/bin/env bash
# Standalone ACMG-AMP classification via fastVEP (Huang-lab/fastVEP).
#
# NOT part of the Nextflow pipeline. Runs after the pipeline's QC step, on the
# QC'd VCFs it publishes (results/<run>/norm_qc/*.norm.vcf), so ACMG is applied
# to exactly the same variants as the ClinVar calls. Emits a per-variant TSV in
# the pipeline's acmg_plp.tsv schema, ready to fold into the carrier matrix later.
#
# Prereqs on Minerva (one-time) — see README.md:
#   * fastvep binary on PATH (cargo install, or conda)
#   * GFF3 gene models + reference FASTA (+ .fai)
#   * fastVEP supplementary DBs built into $SA_DIR (.osa/.oga): gnomAD, ClinVar,
#     REVEL, gene-level constraint/clingen — see fastVEP docs/ACMG_SETUP.md
set -euo pipefail

usage() {
    cat >&2 <<USAGE
Usage: run_fastvep_acmg.sh -i <qc.vcf[.gz]> -o <outdir> \\
         --gff3 <genes.gff3> --fasta <ref.fa> --sa-dir <sa_databases/> \\
         [--acmg-config <config.toml>] [--fastvep <path/to/fastvep>]

Produces:
  <outdir>/<name>.fastvep.vcf        full fastVEP VCF (CSQ incl. ACMG/ACMG_CRITERIA)
  <outdir>/<name>.acmg_plp.tsv       per-variant ACMG table (pipeline schema)
USAGE
    exit 2
}

FASTVEP="${FASTVEP:-fastvep}"
ACMG_CONFIG=""
IN=""; OUTDIR=""; GFF3=""; FASTA=""; SA_DIR=""
while [ $# -gt 0 ]; do
    case "$1" in
        -i|--input) IN="$2"; shift 2 ;;
        -o|--outdir) OUTDIR="$2"; shift 2 ;;
        --gff3) GFF3="$2"; shift 2 ;;
        --fasta) FASTA="$2"; shift 2 ;;
        --sa-dir) SA_DIR="$2"; shift 2 ;;
        --acmg-config) ACMG_CONFIG="$2"; shift 2 ;;
        --fastvep) FASTVEP="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$IN" ] && [ -n "$OUTDIR" ] && [ -n "$GFF3" ] && [ -n "$FASTA" ] && [ -n "$SA_DIR" ] || usage

command -v "$FASTVEP" >/dev/null 2>&1 || { echo "ERROR: fastvep not found ($FASTVEP)" >&2; exit 3; }
for f in "$IN" "$GFF3" "$FASTA"; do [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 3; }; done
[ -d "$SA_DIR" ] || { echo "ERROR: --sa-dir not a directory: $SA_DIR" >&2; exit 3; }

mkdir -p "$OUTDIR"
HERE="$(cd "$(dirname "$0")" && pwd)"
name="$(basename "$IN")"; name="${name%.gz}"; name="${name%.vcf}"; name="${name%.norm}"
VEPOUT="$OUTDIR/${name}.fastvep.vcf"
ACMGOUT="$OUTDIR/${name}.acmg_plp.tsv"

cfg_arg=(); [ -n "$ACMG_CONFIG" ] && cfg_arg=(--acmg-config "$ACMG_CONFIG")

echo "[fastvep] annotating $IN → $VEPOUT" >&2
"$FASTVEP" annotate \
    --input "$IN" \
    --output "$VEPOUT" \
    --gff3 "$GFF3" \
    --fasta "$FASTA" \
    --sa-dir "$SA_DIR" \
    --acmg "${cfg_arg[@]}" \
    --output-format vcf

echo "[parse] $VEPOUT → $ACMGOUT" >&2
python3 "$HERE/parse_fastvep_acmg.py" --vcf "$VEPOUT" --out "$ACMGOUT"

echo "[done] $ACMGOUT" >&2
