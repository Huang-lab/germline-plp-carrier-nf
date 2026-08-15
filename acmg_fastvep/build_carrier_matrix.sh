#!/usr/bin/env bash
# Standalone carrier-matrix assembly (no VEP/fastVEP recompute).
#
# Combines, for ONE chromosome/run, the ClinVar + fastVEP-ACMG per-variant
# classifications with genotypes from the QC'd VCFs into a single per-person
# carrier matrix — reusing the pipeline's bin/build_carrier_matrix.py. This is
# the cheap way to fold standalone fastVEP ACMG results into the ClinVar carrier
# output without re-running the whole Nextflow pipeline.
#
# Inputs (per chromosome) — reuses files you already have:
#   --norm-dir     results-chrN/norm_qc            (*.norm.vcf.gz  → genotypes)
#   --clinvar-dir  results-chrN/variants           (*.clinvar_plp.tsv)
#   --acmg-dir     results-fastvep-chrN            (*.acmg_plp.tsv)
#
# Output: one carrier_matrix.tsv (long): chr,pos,ref,alt,gene,person_id,GT,
#   zygosity,is_clinvar_PLP,is_acmg_PLP,is_AM_PLP  — carriers of any P/LP variant
#   under either framework.
set -euo pipefail

usage() {
    cat >&2 <<'U'
Usage: build_carrier_matrix.sh --norm-dir <norm_qc/> -o <out.tsv> [options]

Required:
  --norm-dir <dir>     dir of QC'd *.norm.vcf.gz (genotype source)
  -o, --out <file>     output combined carrier_matrix.tsv
  (and at least one of:)
  --clinvar-dir <dir>  dir of *.clinvar_plp.tsv   (e.g. results-chrN/variants)
  --acmg-dir <dir>     dir of *.acmg_plp.tsv       (e.g. results-fastvep-chrN)

Options:
  --keep <file>        sample keep-list (one ID per line)
  --bcftools <bin>     bcftools binary (default: bcftools on PATH; else python fallback)
  --out-wide <file>    also write a wide pivot
U
    exit 2
}

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
NORM=""; OUT=""; CVDIR=""; ACDIR=""; KEEP=""; BCF="bcftools"; WIDE=""
while [ $# -gt 0 ]; do
    case "$1" in
        --norm-dir)     NORM="$2"; shift 2 ;;
        -o|--out)       OUT="$2"; shift 2 ;;
        --clinvar-dir)  CVDIR="$2"; shift 2 ;;
        --acmg-dir)     ACDIR="$2"; shift 2 ;;
        --keep)         KEEP="$2"; shift 2 ;;
        --bcftools)     BCF="$2"; shift 2 ;;
        --out-wide)     WIDE="$2"; shift 2 ;;
        -h|--help)      usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$NORM" ] && [ -n "$OUT" ] || usage
[ -d "$NORM" ] || { echo "ERROR: --norm-dir not a directory: $NORM" >&2; exit 3; }
[ -n "$CVDIR$ACDIR" ] || { echo "ERROR: need --clinvar-dir and/or --acmg-dir" >&2; exit 3; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
mkdir -p "$(dirname "$OUT")"

# concat per-chunk TSVs into one, keeping a single header
concat() {   # dir  glob  out
    local dir="$1" g="$2" out="$3" first=1
    : > "$out"
    shopt -s nullglob
    for f in "$dir"/$g; do
        [ -s "$f" ] || continue
        if [ "$first" = 1 ]; then cat "$f" >> "$out"; first=0; else tail -n +2 "$f" >> "$out"; fi
    done
    shopt -u nullglob
    [ -s "$out" ]
}

CV=""; AC=""
[ -n "$CVDIR" ] && concat "$CVDIR" "*.clinvar_plp.tsv" "$WORK/clinvar.tsv" && CV="$WORK/clinvar.tsv" || true
[ -n "$ACDIR" ] && concat "$ACDIR" "*.acmg_plp.tsv"    "$WORK/acmg.tsv"    && AC="$WORK/acmg.tsv"    || true
[ -n "$CV$AC" ] || { echo "ERROR: no classifier TSVs found under the given dirs" >&2; exit 3; }

# qualifying positions = union of is_clinvar_PLP==1 and is_acmg_PLP==1 (header-located columns)
: > "$WORK/qual.txt"
for spec in ${CV:+"$CV:is_clinvar_PLP"} ${AC:+"$AC:is_acmg_PLP"}; do
    f="${spec%%:*}"; flag="${spec##*:}"
    awk -F'\t' -v flag="$flag" '
        NR==1 { for (i=1;i<=NF;i++) c[$i]=i; next }
        (c[flag] && $(c[flag])==1) { print $(c["chr"])"\t"$(c["pos"]) }
    ' "$f" >> "$WORK/qual.txt"
done
sort -u "$WORK/qual.txt" > "$WORK/qual.pos.txt"
echo "[carrier] $(wc -l < "$WORK/qual.pos.txt") qualifying positions" >&2

# genotypes at qualifying positions from every norm VCF → gt.tsv (chr pos ref alt sample GT)
: > "$WORK/gt.tsv"
shopt -s nullglob
vcfs=("$NORM"/*.norm.vcf.gz)
[ ${#vcfs[@]} -eq 0 ] && vcfs=("$NORM"/*.vcf.gz)
[ ${#vcfs[@]} -eq 0 ] && vcfs=("$NORM"/*.norm.vcf)
[ ${#vcfs[@]} -eq 0 ] && vcfs=("$NORM"/*.vcf)
shopt -u nullglob
[ ${#vcfs[@]} -gt 0 ] || { echo "ERROR: no VCFs in $NORM" >&2; exit 3; }

if [ -s "$WORK/qual.pos.txt" ]; then
    have_bcf=0; command -v "$BCF" >/dev/null 2>&1 && have_bcf=1
    for vcf in "${vcfs[@]}"; do
        [ -s "$vcf" ] || continue
        if [ "$have_bcf" = 1 ]; then
            "$BCF" query -T "$WORK/qual.pos.txt" \
                -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' "$vcf" \
              | awk -F'\t' 'BEGIN{OFS="\t"}{for(i=5;i<=NF;i++){split($i,kv,"=");print $1,$2,$3,$4,kv[1],kv[2]}}' \
              >> "$WORK/gt.tsv"
        else
            python3 "$REPO/bin/fixture_gt_extract.py" --vcf "$vcf" --positions "$WORK/qual.pos.txt" --out "$WORK/gt.chunk"
            cat "$WORK/gt.chunk" >> "$WORK/gt.tsv"
        fi
    done
fi

python3 "$REPO/bin/build_carrier_matrix.py" \
    --gt "$WORK/gt.tsv" \
    ${CV:+--clinvar "$CV"} ${AC:+--acmg "$AC"} \
    ${KEEP:+--keep "$KEEP"} \
    --out-long "$OUT" ${WIDE:+--out-wide "$WIDE"}

echo "[done] $OUT  ($(( $(wc -l < "$OUT") - 1 )) carrier rows)" >&2
