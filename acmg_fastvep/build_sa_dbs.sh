#!/usr/bin/env bash
# Build fastVEP supplementary-annotation databases into $FV_ROOT/sa_db.
#
# Run on a COMPUTE node (offline — reads local source files in $FV_ROOT/sa_src).
# Whole-genome sources (clinvar, revel, gnomad_genes) are built once if missing.
# gnomAD is per-chromosome; pass chromosomes to build, or none to auto-detect
# from the gnomad.exomes.*.sites.chrN.*.bgz files present in sa_src.
#
# Usage:
#   export FASTVEP=/.../fastvep/bin/fastvep
#   export FV_ROOT=/.../data/fastvep            # has sa_src/ ; sa_db/ is created
#   export CLINVAR=/.../clinvar_YYYYMMDD.vcf.gz  # optional (whole-genome clinvar source)
#   acmg_fastvep/build_sa_dbs.sh                 # build whole-genome + all gnomAD chrs found
#   acmg_fastvep/build_sa_dbs.sh 2 3 4           # build whole-genome (if missing) + gnomAD chr2,3,4
set -euo pipefail

FASTVEP="${FASTVEP:?set FASTVEP=<fastvep binary>}"
FV_ROOT="${FV_ROOT:?set FV_ROOT=<tool root with sa_src/>}"
SRC="$FV_ROOT/sa_src"
DB="$FV_ROOT/sa_db"
CLINVAR="${CLINVAR:-}"
mkdir -p "$DB"

command -v "$FASTVEP" >/dev/null 2>&1 || { echo "ERROR: fastvep not found: $FASTVEP" >&2; exit 3; }
[ -d "$SRC" ] || { echo "ERROR: sources dir missing: $SRC" >&2; exit 3; }

build_once() {   # source  outbase  inputfile
    local src="$1" out="$2" in="$3"
    if ls "$DB/$out".* >/dev/null 2>&1; then echo "[skip] $out already built"; return; fi
    [ -s "$in" ] || { echo "[skip] $out — no input file: $in" >&2; return; }
    echo "[build] $out  <-  $in"
    "$FASTVEP" sa-build --source "$src" -i "$in" -o "$DB/$out" --assembly GRCh38
}

# --- whole-genome sources (built once) ---
[ -n "$CLINVAR" ] && build_once clinvar clinvar "$CLINVAR"
build_once revel revel "$SRC/revel_with_transcript_ids"
GC="$(ls "$SRC"/gnomad*constraint*metrics*.tsv 2>/dev/null | head -1 || true)"
[ -n "$GC" ] && build_once gnomad_genes gnomad_genes "$GC"

# --- gnomAD per-chromosome ---
chrs=("$@")
if [ ${#chrs[@]} -eq 0 ]; then
    shopt -s nullglob
    for f in "$SRC"/gnomad.exomes.*.sites.chr*.vcf.bgz; do
        c="$(basename "$f" | sed -E 's/.*\.sites\.chr([0-9XY]+)\.vcf\.bgz/\1/')"
        chrs+=("$c")
    done
    shopt -u nullglob
fi
for c in "${chrs[@]}"; do
    f="$(ls "$SRC"/gnomad.exomes.*.sites.chr${c}.vcf.bgz 2>/dev/null | head -1 || true)"
    [ -n "$f" ] || { echo "[skip] no gnomAD source for chr$c in $SRC" >&2; continue; }
    build_once gnomad "gnomad_chr${c}" "$f"
done

echo "[done] databases in $DB:"
ls -lh "$DB"
