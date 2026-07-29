#!/usr/bin/env bash
# Standalone fastVEP runner for a SINGLE QC'd VCF (Huang-lab/fastVEP).
#
# NOT part of the Nextflow pipeline. Runs after the pipeline's QC step, on the
# QC'd VCFs it publishes (results/<run>/norm_qc/<chunk>.norm.vcf.gz), so fastVEP
# annotates exactly the same variants as the ClinVar path — without touching the
# working Nextflow pipeline.
#
# Two modes:
#   * plain annotation (default)  — needs only a GFF3 gene model (+ optional
#     FASTA for HGVS). No supplementary databases. This is the "just make it run"
#     path. Emits <name>.fastvep.vcf.
#   * ACMG (--acmg)               — additionally requires --sa-dir pointing at
#     fastVEP supplementary DBs (.osa/.oga: gnomAD, ClinVar, REVEL, gene-level;
#     see fastVEP docs/ACMG_SETUP.md). Emits <name>.fastvep.vcf AND
#     <name>.acmg_plp.tsv (pipeline acmg_plp.tsv schema) via parse_fastvep_acmg.py.
#
# Prereqs: `fastvep` on PATH (see setup_fastvep.sh). fastVEP is a native binary —
# no container / Singularity module needed (unlike the Nextflow steps).
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: run_fastvep.sh -i <qc.vcf[.gz]> -o <outdir> --gff3 <genes.gff3> [options]

Required:
  -i, --input <vcf>       QC'd VCF (.vcf or .vcf.gz), Ensembl contigs
  -o, --outdir <dir>      output directory (created if absent)
  --gff3 <gff3>           Ensembl GFF3 gene model

Options:
  --fasta <fa>            reference FASTA (needed for --hgvs / sequence context)
  --hgvs                  emit HGVSc/HGVSp (requires --fasta)
  --pick                  most-severe consequence per variant (smaller output)
  --acmg                  run ACMG-AMP classification (REQUIRES --sa-dir)
  --sa-dir <dir>          fastVEP supplementary-annotation DB directory
  --acmg-config <toml>    ACMG threshold overrides
  --threads <N>           worker threads (exported as RAYON_NUM_THREADS)
  --fastvep <path>        path to the fastvep binary (default: fastvep, or $FASTVEP)
  -h, --help

Outputs:
  <outdir>/<name>.fastvep.vcf     annotated VCF (CSQ; ACMG cols filled only with --acmg)
  <outdir>/<name>.acmg_plp.tsv    per-variant ACMG table (only in --acmg mode)
USAGE
    exit 2
}

FASTVEP="${FASTVEP:-fastvep}"
IN=""; OUTDIR=""; GFF3=""; FASTA=""; SA_DIR=""; ACMG_CONFIG=""; THREADS=""
DO_ACMG=0; DO_HGVS=0; DO_PICK=0
while [ $# -gt 0 ]; do
    case "$1" in
        -i|--input)     IN="$2"; shift 2 ;;
        -o|--outdir)    OUTDIR="$2"; shift 2 ;;
        --gff3)         GFF3="$2"; shift 2 ;;
        --fasta)        FASTA="$2"; shift 2 ;;
        --sa-dir)       SA_DIR="$2"; shift 2 ;;
        --acmg-config)  ACMG_CONFIG="$2"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        --fastvep)      FASTVEP="$2"; shift 2 ;;
        --acmg)         DO_ACMG=1; shift ;;
        --hgvs)         DO_HGVS=1; shift ;;
        --pick)         DO_PICK=1; shift ;;
        -h|--help)      usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done

# --- validation ---
[ -n "$IN" ] && [ -n "$OUTDIR" ] && [ -n "$GFF3" ] || usage
command -v "$FASTVEP" >/dev/null 2>&1 || { echo "ERROR: fastvep not found ($FASTVEP) — see setup_fastvep.sh" >&2; exit 3; }
for f in "$IN" "$GFF3"; do [ -s "$f" ] || { echo "ERROR: missing/empty: $f" >&2; exit 3; }; done
if [ -n "$FASTA" ]; then [ -s "$FASTA" ] || { echo "ERROR: missing/empty FASTA: $FASTA" >&2; exit 3; }; fi
if [ "$DO_HGVS" = 1 ] && [ -z "$FASTA" ]; then echo "ERROR: --hgvs requires --fasta" >&2; exit 3; fi
if [ "$DO_ACMG" = 1 ]; then
    [ -n "$SA_DIR" ] || { echo "ERROR: --acmg requires --sa-dir (supplementary DBs)" >&2; exit 3; }
    [ -d "$SA_DIR" ] || { echo "ERROR: --sa-dir not a directory: $SA_DIR" >&2; exit 3; }
fi
[ -n "$THREADS" ] && export RAYON_NUM_THREADS="$THREADS"

mkdir -p "$OUTDIR"
HERE="$(cd "$(dirname "$0")" && pwd)"
name="$(basename "$IN")"; name="${name%.gz}"; name="${name%.vcf}"; name="${name%.norm}"
VEPOUT="$OUTDIR/${name}.fastvep.vcf"

# --- build the annotate command ---
cmd=("$FASTVEP" annotate --input "$IN" --output "$VEPOUT" --gff3 "$GFF3" --output-format vcf)
[ -n "$FASTA" ]      && cmd+=(--fasta "$FASTA")
[ "$DO_HGVS" = 1 ]   && cmd+=(--hgvs)
[ "$DO_PICK" = 1 ]   && cmd+=(--pick)
if [ "$DO_ACMG" = 1 ]; then
    cmd+=(--sa-dir "$SA_DIR" --acmg)
    [ -n "$ACMG_CONFIG" ] && cmd+=(--acmg-config "$ACMG_CONFIG")
fi

echo "[fastvep] ${cmd[*]}" >&2
"${cmd[@]}"

# --- ACMG post-parse (only in --acmg mode) ---
if [ "$DO_ACMG" = 1 ]; then
    ACMGOUT="$OUTDIR/${name}.acmg_plp.tsv"
    echo "[parse] $VEPOUT -> $ACMGOUT" >&2
    python3 "$HERE/parse_fastvep_acmg.py" --vcf "$VEPOUT" --out "$ACMGOUT"
    echo "[done] $VEPOUT  $ACMGOUT" >&2
else
    echo "[done] $VEPOUT" >&2
fi
