#!/usr/bin/env bash
# Batch driver: run fastVEP over EVERY QC'd VCF in a norm_qc/ directory.
#
# For each <chunk>.norm.vcf.gz it calls run_fastvep.sh (plain or --acmg). Runs
# locally in sequence by default, or submits one LSF job per chunk with --lsf.
# Writes a chunk->status manifest so a failed chunk is never silent (the
# lesson from the batch3 dropped-chunk incident). In --acmg mode, optionally
# concatenates the per-chunk *.acmg_plp.tsv into one table.
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: run_fastvep_batch.sh --in-dir <norm_qc/> -o <outdir> --gff3 <g> [options]

Required:
  --in-dir <dir>          directory of QC'd VCFs (*.norm.vcf.gz / *.vcf.gz)
  -o, --outdir <dir>      output directory
  --gff3 <gff3>           Ensembl GFF3

Passthrough to run_fastvep.sh:
  --fasta <fa>  --hgvs  --pick  --acmg  --sa-dir <dir>  --acmg-config <toml>
  --threads <N>  --fastvep <path>

Batch options:
  --lsf                   submit one bsub per chunk (else run locally, sequential)
  -P, --project <alloc>   LSF allocation (default: $MINERVA_ALLOCATION)
  --queue <q>             LSF queue (default: premium)
  --walltime <hh:mm>      LSF walltime (default: 4:00)
  --concat                after (local) runs, concatenate *.acmg_plp.tsv into one
  -h, --help
USAGE
    exit 2
}

HERE="$(cd "$(dirname "$0")" && pwd)"
INDIR=""; OUTDIR=""; GFF3=""; FASTA=""; SA_DIR=""; ACMG_CONFIG=""; THREADS=""; FASTVEP=""
DO_ACMG=0; DO_HGVS=0; DO_PICK=0; USE_LSF=0; DO_CONCAT=0
PROJECT="${MINERVA_ALLOCATION:-}"; QUEUE="premium"; WALL="4:00"
while [ $# -gt 0 ]; do
    case "$1" in
        --in-dir)       INDIR="$2"; shift 2 ;;
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
        --lsf)          USE_LSF=1; shift ;;
        -P|--project)   PROJECT="$2"; shift 2 ;;
        --queue)        QUEUE="$2"; shift 2 ;;
        --walltime)     WALL="$2"; shift 2 ;;
        --concat)       DO_CONCAT=1; shift ;;
        -h|--help)      usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$INDIR" ] && [ -n "$OUTDIR" ] && [ -n "$GFF3" ] || usage
[ -d "$INDIR" ] || { echo "ERROR: --in-dir not a directory: $INDIR" >&2; exit 3; }
if [ "$USE_LSF" = 1 ] && [ -z "$PROJECT" ]; then
    echo "ERROR: --lsf needs -P/--project or \$MINERVA_ALLOCATION" >&2; exit 3
fi
mkdir -p "$OUTDIR"

# common run_fastvep.sh args shared by every chunk
common=(--gff3 "$GFF3")
[ -n "$FASTA" ]        && common+=(--fasta "$FASTA")
[ "$DO_HGVS" = 1 ]     && common+=(--hgvs)
[ "$DO_PICK" = 1 ]     && common+=(--pick)
[ "$DO_ACMG" = 1 ]     && common+=(--acmg)
[ -n "$SA_DIR" ]       && common+=(--sa-dir "$SA_DIR")
[ -n "$ACMG_CONFIG" ]  && common+=(--acmg-config "$ACMG_CONFIG")
[ -n "$THREADS" ]      && common+=(--threads "$THREADS")
[ -n "$FASTVEP" ]      && common+=(--fastvep "$FASTVEP")

# collect inputs (norm.vcf.gz preferred; fall back to any vcf.gz / vcf)
shopt -s nullglob
inputs=("$INDIR"/*.norm.vcf.gz)
[ ${#inputs[@]} -eq 0 ] && inputs=("$INDIR"/*.vcf.gz)
[ ${#inputs[@]} -eq 0 ] && inputs=("$INDIR"/*.vcf)
shopt -u nullglob
[ ${#inputs[@]} -gt 0 ] || { echo "ERROR: no VCFs found in $INDIR" >&2; exit 3; }
echo "[batch] ${#inputs[@]} input VCF(s) from $INDIR" >&2

MANIFEST="$OUTDIR/batch_manifest.tsv"
: > "$MANIFEST"
n_submitted=0; n_ok=0; n_fail=0; n_skipped=0
for vcf in "${inputs[@]}"; do
    name="$(basename "$vcf")"; name="${name%.gz}"; name="${name%.vcf}"; name="${name%.norm}"

    # Skip a chunk whose expected output(s) already exist and are non-empty —
    # makes the batch driver safe to re-run after a partial/interrupted run
    # (e.g. resubmitting after adding new SA databases, without redoing chunks
    # that already completed under the old ones).
    vep_out="$OUTDIR/${name}.fastvep.vcf.gz"
    acmg_out="$OUTDIR/${name}.acmg_plp.tsv"
    already_done=0
    if [ -s "$vep_out" ]; then
        if [ "$DO_ACMG" = 1 ]; then
            [ -s "$acmg_out" ] && already_done=1
        else
            already_done=1
        fi
    fi
    if [ "$already_done" = 1 ]; then
        printf '%s\tskipped-already-done\n' "$name" >> "$MANIFEST"
        n_skipped=$((n_skipped+1))
        continue
    fi

    if [ "$USE_LSF" = 1 ]; then
        bsub -P "$PROJECT" -q "$QUEUE" -n "${THREADS:-2}" \
             -R "rusage[mem=8000]" -W "$WALL" \
             -o "$OUTDIR/${name}.fastvep.%J.out" -e "$OUTDIR/${name}.fastvep.%J.err" \
             "$HERE/run_fastvep.sh" -i "$vcf" -o "$OUTDIR" "${common[@]}"
        printf '%s\tsubmitted\n' "$name" >> "$MANIFEST"
        n_submitted=$((n_submitted+1))
    else
        if "$HERE/run_fastvep.sh" -i "$vcf" -o "$OUTDIR" "${common[@]}"; then
            printf '%s\tok\n' "$name" >> "$MANIFEST"; n_ok=$((n_ok+1))
        else
            printf '%s\tFAILED\n' "$name" >> "$MANIFEST"; n_fail=$((n_fail+1))
            echo "[batch] WARN: chunk FAILED: $name" >&2
        fi
    fi
done

if [ "$USE_LSF" = 1 ]; then
    echo "[batch] submitted $n_submitted LSF jobs, skipped $n_skipped already-done. Manifest: $MANIFEST" >&2
    echo "[batch] (concat, if wanted, after jobs finish: rerun with --concat and no --lsf, or concat manually)" >&2
    exit 0
fi

echo "[batch] done: $n_ok ok, $n_fail failed, $n_skipped skipped (already done) of ${#inputs[@]}. Manifest: $MANIFEST" >&2

# optional concat of per-chunk ACMG tables into one
if [ "$DO_CONCAT" = 1 ] && [ "$DO_ACMG" = 1 ]; then
    OUT_ALL="$OUTDIR/all.acmg_plp.tsv"
    first=1
    : > "$OUT_ALL"
    for tsv in "$OUTDIR"/*.acmg_plp.tsv; do
        [ "$tsv" = "$OUT_ALL" ] && continue
        [ -s "$tsv" ] || continue
        if [ "$first" = 1 ]; then cat "$tsv" >> "$OUT_ALL"; first=0
        else tail -n +2 "$tsv" >> "$OUT_ALL"; fi
    done
    echo "[batch] concatenated ACMG table: $OUT_ALL" >&2
fi

[ "$n_fail" -eq 0 ] || exit 1
