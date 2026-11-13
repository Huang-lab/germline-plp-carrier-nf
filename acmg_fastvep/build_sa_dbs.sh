#!/usr/bin/env bash
# Build the fastVEP v0.5.0 ACMG supplementary-annotation stack from
# $FV_HOME/sa_src (fetch_sa_sources.sh) into $FV_SA (= $FV_HOME/sa_db), exactly
# as fastVEP docs/ACMG_SETUP.md at the v0.5.0 tag builds it. Offline: run it on
# COMPUTE nodes, via build_sa_dbs.lsf.
#
# Two modes, so the work splits into one job plus one array:
#   --single     the five whole-genome databases
#                  clinvar.osa2          <- clinvar_<date>.vcf.gz
#                  clinvar_protein.oga   <- variant_summary.txt.gz (NOT the VCF: no splice index)
#                  gnomad_genes.oga      <- gnomad.v4.1.constraint_metrics.tsv
#                  omim.oga              <- ClinGen GDV via clingen_gdv_to_oga.py
#                                           (real OMIM genemap2.txt only with --use-omim)
#                  repeatmasker.osi      <- rmsk.txt.gz via repeatmasker_to_bed.py,
#                                           --source custom_bed --name repeatmasker
#   --chrom N    the four per-chromosome databases for chromosome N (1-22, X, Y)
#                  gnomad_chrN.osa2      <- gnomAD v4.1 exomes chrN, --chunk-bits 16
#                  spliceai_chrN.osa2    <- extract_gnomad_scores.py (spliceai_ds_max), --chunk-bits 16
#                  phylop_chrN.osa2      <- extract_gnomad_scores.py (phylop),          --chunk-bits 16
#                  revel_chrN.osa2       <- sa_src/revel/revel_chrN.csv (none for Y: REVEL has no chrY)
#
# Extensions are what v0.5.0's `sa-build --format auto` (the default) writes:
# .osa2 for every allele-level source (OSA2_SUPPORTED_SOURCES in
# crates/fastvep-cli/src/pipeline/sa_build.rs includes clinvar and revel), .oga
# for the gene-level sources, .osi for custom_bed. (ACMG_SETUP.md's size table
# lists clinvar.osa / revel_chr*.osa from the benchmark's older build; auto
# writes .osa2.)
#
# Safety rules, each one a past or documented silent failure:
#   * Never --allow-empty: an empty database answers every query "no record",
#     which PM2 reads as "absent from the population".
#   * A database is built in $FV_HOME/sa_build_tmp and moved into sa_db only
#     after sa-build exits 0, then <db>.done is written. A killed job leaves
#     nothing in sa_db for fastVEP to half-load.
#   * Skip = output AND .done both present. Either alone is rebuilt.
#   * Only gnomAD v4.1 EXOMES are read. Never put joint or genomes databases
#     beside these: fastVEP keeps one gnomad answer per allele (the database
#     sorting last), so mixing silently swaps frequencies.
#   * The binary must be the pinned one (sha256 == FV_BIN_SHA256 from env.sh).
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: build_sa_dbs.sh (--single | --chrom <N>) [options]

Environment (from `source $FV_HOME/env.sh`, or FV_ENV=<path to env.sh>):
  FV_HOME FV_BIN FV_BIN_SHA256 FV_SRC [FV_SA, default $FV_HOME/sa_db]

Options:
  --single              build clinvar, clinvar_protein, gnomad_genes, omim, repeatmasker
  --chrom <N>           build gnomad/spliceai/phylop/revel for chromosome N (1-22, X, Y)
  --use-omim            omim.oga from sa_src/genemap2.txt (licensed OMIM) instead of ClinGen
  --delete-gnomad-src   after gnomad+spliceai+phylop for --chrom N all succeed, delete
                        that chromosome's gnomAD .vcf.bgz/.tbi from sa_src
  --force               rebuild even when <db> and <db>.done exist
  -h, --help
USAGE
    exit 2
}

if [ -n "${FV_ENV:-}" ]; then
    # shellcheck disable=SC1090
    . "$FV_ENV"
fi
MODE=""; CHROM=""; USE_OMIM=0; DELETE_GNOMAD=0; FORCE=0
while [ $# -gt 0 ]; do
    case "$1" in
        --single)            MODE=single; shift ;;
        --chrom)             MODE=chrom; CHROM="${2#chr}"; shift 2 ;;
        --use-omim)          USE_OMIM=1; shift ;;
        --delete-gnomad-src) DELETE_GNOMAD=1; shift ;;
        --force)             FORCE=1; shift ;;
        -h|--help)           usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$MODE" ] || usage
if [ "$MODE" = chrom ]; then
    case "$CHROM" in
        [1-9]|1[0-9]|2[0-2]|X|Y) ;;
        *) echo "ERROR: --chrom must be 1-22, X or Y (got '$CHROM')" >&2; exit 2 ;;
    esac
fi

: "${FV_HOME:?set FV_HOME (source \$FV_HOME/env.sh, or FV_ENV=<env.sh>)}"
: "${FV_BIN:?set FV_BIN (source \$FV_HOME/env.sh)}"
: "${FV_BIN_SHA256:?set FV_BIN_SHA256 (source \$FV_HOME/env.sh)}"
: "${FV_SRC:?set FV_SRC (source \$FV_HOME/env.sh)}"
SRC="$FV_HOME/sa_src"
DB="${FV_SA:-$FV_HOME/sa_db}"
TMPROOT="$FV_HOME/sa_build_tmp"
LOGDIR="$FV_HOME/sa_build_logs"
mkdir -p "$DB" "$TMPROOT" "$LOGDIR" "$SRC/derived"

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}
die() { echo "ERROR: $*" >&2; exit 4; }

[ -x "$FV_BIN" ] || die "FV_BIN not executable: $FV_BIN"
[ "$(sha256_of "$FV_BIN")" = "$FV_BIN_SHA256" ] || \
    die "sha256($FV_BIN) != FV_BIN_SHA256 from env.sh — not the pinned binary; re-run setup_fastvep.sh"
[ -d "$SRC" ] || die "sources dir missing: $SRC (run fetch_sa_sources.sh)"
SCRIPTS="$FV_SRC/scripts"
SA_SOURCES="$FV_SRC/analysis/acmg_benchmark/scripts/sa_sources"

WORK="$(mktemp -d "$TMPROOT/${MODE}${CHROM}.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

BUILT=""; SKIPPED=""
# build <source> <db base> <ext> <input> [extra sa-build args...]
build() {
    local source="$1" base="$2" ext="$3" in="$4"; shift 4
    local out="$DB/$base.$ext" done="$DB/$base.$ext.done" log="$LOGDIR/$base.log"
    if [ "$FORCE" = 0 ] && [ -s "$out" ] && [ -e "$done" ]; then
        echo "[skip] $base.$ext (present with .done)" >&2
        SKIPPED="$SKIPPED $base"; return 0
    fi
    [ -s "$in" ] || die "$base: input missing or empty: $in"
    rm -f "$out" "$done"
    mkdir -p "$WORK/$base"
    echo "[build] $base.$ext  <-  $in" >&2
    if ! "$FV_BIN" sa-build --source "$source" -i "$in" -o "$WORK/$base/$base" \
            --assembly GRCh38 --no-progress "$@" > "$log" 2>&1; then
        cat "$log" >&2
        die "$base: fastvep sa-build failed (log: $log); nothing published"
    fi
    [ -s "$WORK/$base/$base.$ext" ] || {
        cat "$log" >&2
        ls -l "$WORK/$base" >&2
        die "$base: sa-build exited 0 but wrote no $base.$ext (extension changed upstream?)"
    }
    mv -f "$WORK/$base/$base.$ext" "$out"
    local in_sha=""
    [ -s "$in.sha256" ] && in_sha="$(cat "$in.sha256")"
    {
        echo "built_utc     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "source        $source"
        echo "input         $in"
        echo "input_sha256  ${in_sha:-not recorded (derived file)}"
        echo "args          --assembly GRCh38 --no-progress $*"
        echo "fastvep       $FV_BIN"
        echo "fastvep_sha256 $FV_BIN_SHA256"
        echo "output_sha256 $(sha256_of "$out")"
        echo "host          $(hostname 2>/dev/null || echo unknown)  LSB_JOBID=${LSB_JOBID:-}"
    } > "$done"
    rm -rf "${WORK:?}/$base"
    BUILT="$BUILT $base"
}
is_done() { [ -s "$DB/$1" ] && [ -e "$DB/$1.done" ]; }

if [ "$MODE" = single ]; then
    # clinvar: the dated release fetch_sa_sources.sh recorded and md5-verified.
    [ -s "$SRC/CLINVAR_DATE" ] || die "no $SRC/CLINVAR_DATE (run fetch_sa_sources.sh)"
    CV="$SRC/clinvar_$(cat "$SRC/CLINVAR_DATE").vcf.gz"
    [ -e "$CV.md5.ok" ] || die "$CV was not md5-verified by fetch_sa_sources.sh"
    build clinvar clinvar osa2 "$CV"

    build clinvar_protein clinvar_protein oga "$SRC/variant_summary.txt.gz"
    build gnomad_genes gnomad_genes oga "$SRC/gnomad.v4.1.constraint_metrics.tsv"

    # omim (json key) = ClinGen Gene-Disease Validity by default, per ACMG_SETUP.md.
    if [ "$USE_OMIM" = 1 ]; then
        build omim omim oga "$SRC/genemap2.txt"
    elif is_done omim.oga && [ "$FORCE" = 0 ]; then
        build omim omim oga "$SRC/derived/clingen_gdv.tsv"   # reports the skip
    else
        shopt -s nullglob; cg=("$SRC"/clingen_gene_validity_*.csv); shopt -u nullglob
        [ ${#cg[@]} -eq 1 ] || die "expected exactly one clingen_gene_validity_*.csv in $SRC, found ${#cg[@]}"
        # Two positional arguments: input CSV, output TSV.
        python3 "$SA_SOURCES/clingen_gdv_to_oga.py" "${cg[0]}" "$WORK/clingen_gdv.tsv" \
            2> "$LOGDIR/clingen_gdv_to_oga.log" || { cat "$LOGDIR/clingen_gdv_to_oga.log" >&2; die "clingen_gdv_to_oga.py failed"; }
        cat "$LOGDIR/clingen_gdv_to_oga.log" >&2
        cp -f "$WORK/clingen_gdv.tsv" "$SRC/derived/clingen_gdv.tsv"   # check_sa_db.sh reads it
        build omim omim oga "$SRC/derived/clingen_gdv.tsv"
    fi

    if is_done repeatmasker.osi && [ "$FORCE" = 0 ]; then
        echo "[skip] repeatmasker.osi (present with .done)" >&2; SKIPPED="$SKIPPED repeatmasker"
    else
        [ -s "$SRC/rmsk.txt.gz" ] || die "repeatmasker: missing $SRC/rmsk.txt.gz"
        python3 "$SA_SOURCES/repeatmasker_to_bed.py" "$SRC/rmsk.txt.gz" > "$WORK/repeatmasker.bed" \
            2> "$LOGDIR/repeatmasker_to_bed.log" || { cat "$LOGDIR/repeatmasker_to_bed.log" >&2; die "repeatmasker_to_bed.py failed"; }
        cat "$LOGDIR/repeatmasker_to_bed.log" >&2
        # --name repeatmasker is load-bearing: the classifier finds BP3's track by a key containing "repeat".
        build custom_bed repeatmasker osi "$WORK/repeatmasker.bed" --name repeatmasker
    fi
else
    c="$CHROM"
    G="$SRC/gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz"
    build gnomad "gnomad_chr${c}" osa2 "$G" --chunk-bits 16

    # SpliceAI and PhyloP are distilled from the same gnomAD file (its
    # spliceai_ds_max and phylop INFO fields). --url points the extractor at the
    # local .bgz; with no --regions it reads the whole chromosome via tabix.
    if [ "$FORCE" = 0 ] && is_done "spliceai_chr${c}.osa2" && is_done "phylop_chr${c}.osa2"; then
        echo "[skip] spliceai_chr${c}.osa2 phylop_chr${c}.osa2 (present with .done)" >&2
        SKIPPED="$SKIPPED spliceai_chr${c} phylop_chr${c}"
    else
        [ -s "$G" ] && [ -s "$G.tbi" ] || die "chr$c: need $G and its .tbi for extract_gnomad_scores.py"
        echo "[extract] spliceai/phylop chr$c <- $G" >&2
        python3 "$SCRIPTS/extract_gnomad_scores.py" --chrom "$c" --url "$G" --out-dir "$WORK/extract" \
            2> "$LOGDIR/extract_gnomad_scores_chr${c}.log" || {
            cat "$LOGDIR/extract_gnomad_scores_chr${c}.log" >&2; die "extract_gnomad_scores.py chr$c failed"; }
        cat "$LOGDIR/extract_gnomad_scores_chr${c}.log" >&2
        build spliceai "spliceai_chr${c}" osa2 "$WORK/extract/spliceai_chr${c}.vcf" --chunk-bits 16
        build phylop "phylop_chr${c}" osa2 "$WORK/extract/phylop_chr${c}.tsv" --chunk-bits 16
    fi

    if [ "$c" = Y ]; then
        echo "[skip] revel_chrY: REVEL v1.3 has no chrY records, and sa-build refuses an empty database" >&2
    else
        build revel "revel_chr${c}" osa2 "$SRC/revel/revel_chr${c}.csv"
    fi

    if [ "$DELETE_GNOMAD" = 1 ]; then
        if is_done "gnomad_chr${c}.osa2" && is_done "spliceai_chr${c}.osa2" && is_done "phylop_chr${c}.osa2"; then
            echo "[clean] --delete-gnomad-src: removing $G(.tbi)" >&2
            rm -f "$G" "$G.tbi"
        else
            echo "[clean] NOT deleting $G: not every database built from it has a .done" >&2
        fi
    fi
fi

echo "[done] mode=$MODE${CHROM:+ chr$CHROM} built:${BUILT:- none} skipped:${SKIPPED:- none}" >&2
echo "[done] verify the whole stack with: acmg_fastvep/check_sa_db.sh" >&2
