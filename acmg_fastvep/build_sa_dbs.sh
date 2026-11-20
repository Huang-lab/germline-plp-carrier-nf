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
#   * Skip = output AND .done present AND the .done records the same upstream
#     source (sha256), binary and sa-build args. Anything else is rebuilt, with
#     the reason printed (a new ClinGen export, --use-omim, a new binary ...).
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

# Identity of an upstream source file: its sha256, cached in <file>.sha256
# with the size and mtime it was taken at in <file>.sha256.stat (the same cache
# fetch_sa_sources.sh keeps); a file whose size or mtime moved is re-hashed.
# The cache also outlives a gnomAD file removed by --delete-gnomad-src, so its
# databases still count as current.
stat_of() { printf '%s %s' "$(wc -c < "$1" | tr -d ' ')" "$(date -r "$1" +%s)"; }
ident() {   # file
    local f="$1"
    if [ ! -e "$f" ] && [ -s "$f.sha256" ]; then cat "$f.sha256"; return 0; fi
    [ -s "$f" ] || die "source missing or empty: $f"
    if [ -s "$f.sha256" ] && [ "$(cat "$f.sha256.stat" 2>/dev/null)" = "$(stat_of "$f")" ]; then
        cat "$f.sha256"; return 0
    fi
    sha256_of "$f" > "$f.sha256"
    stat_of "$f" > "$f.sha256.stat"
    cat "$f.sha256"
}
kv() { awk -v k="$1" '$1 == k && !done { sub(/^[^ ]+ +/, ""); print; done = 1 }' "$2"; }

# up_to_date <db file name> <source file> [sa-build args...]
# The database counts as built only if it and its .done exist AND the .done
# records this source (by sha256), this binary and these sa-build args. A
# changed ClinGen export, a new REVEL split or --use-omim therefore rebuilds
# without --force; the reason is printed.
up_to_date() {
    local db="$1" src="$2"; shift 2
    local done="$DB/$db.done" want_args="--assembly GRCh38 --no-progress $*" why=""
    [ "$FORCE" = 0 ] || return 1
    [ -s "$DB/$db" ] && [ -e "$done" ] || return 1
    [ "$(kv source_file "$done")" = "$src" ] || why="source file was $(kv source_file "$done")"
    [ -n "$why" ] || [ "$(kv source_sha256 "$done")" = "$(ident "$src")" ] || why="source changed (sha256)"
    [ -n "$why" ] || [ "$(kv fastvep_sha256 "$done")" = "$FV_BIN_SHA256" ] || why="built by another fastvep"
    [ -n "$why" ] || [ "$(kv args "$done")" = "$want_args" ] || why="sa-build args differ"
    if [ -n "$why" ]; then
        echo "[rebuild] $db: $why" >&2
        return 1
    fi
    return 0
}

# build <source> <db base> <ext> <input> <upstream source file> [extra sa-build args...]
# <input> is what sa-build reads (it may be derived, e.g. the ClinGen TSV);
# <upstream source file> is what sa_src holds and what the .done is keyed on.
build() {
    local source="$1" base="$2" ext="$3" in="$4" upstream="$5"; shift 5
    local out="$DB/$base.$ext" done="$DB/$base.$ext.done" log="$LOGDIR/$base.log"
    if up_to_date "$base.$ext" "$upstream" "$@"; then
        echo "[skip] $base.$ext (current: same source, binary and args)" >&2
        SKIPPED="$SKIPPED $base"; return 0
    fi
    [ -s "$in" ] || die "$base: input missing or empty: $in"
    local src_sha
    src_sha="$(ident "$upstream")"
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
    {
        echo "built_utc     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "source        $source"
        echo "input         $in"
        echo "source_file   $upstream"
        echo "source_sha256 $src_sha"
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
    build clinvar clinvar osa2 "$CV" "$CV"

    # variant_summary is NCBI's moving file; the copy in sa_src is the pinned one
    # (its date: sa_src/VARIANT_SUMMARY_DATE), and the .done keys on its sha256.
    VS="$SRC/variant_summary.txt.gz"
    build clinvar_protein clinvar_protein oga "$VS" "$VS"
    GC="$SRC/gnomad.v4.1.constraint_metrics.tsv"
    build gnomad_genes gnomad_genes oga "$GC" "$GC"

    # omim (json key) = ClinGen Gene-Disease Validity by default, per ACMG_SETUP.md;
    # --use-omim keys it on genemap2.txt instead, so switching rebuilds by itself.
    if [ "$USE_OMIM" = 1 ]; then
        build omim omim oga "$SRC/genemap2.txt" "$SRC/genemap2.txt"
    else
        shopt -s nullglob; cg=("$SRC"/clingen_gene_validity_*.csv); shopt -u nullglob
        [ ${#cg[@]} -eq 1 ] || die "expected exactly one clingen_gene_validity_*.csv in $SRC, found ${#cg[@]}"
        if up_to_date omim.oga "${cg[0]}"; then
            echo "[skip] omim.oga (current: same ClinGen export, binary and args)" >&2
            SKIPPED="$SKIPPED omim"
        else
            # Two positional arguments: input CSV, output TSV.
            python3 "$SA_SOURCES/clingen_gdv_to_oga.py" "${cg[0]}" "$WORK/clingen_gdv.tsv" \
                2> "$LOGDIR/clingen_gdv_to_oga.log" || { cat "$LOGDIR/clingen_gdv_to_oga.log" >&2; die "clingen_gdv_to_oga.py failed"; }
            cat "$LOGDIR/clingen_gdv_to_oga.log" >&2
            cp -f "$WORK/clingen_gdv.tsv" "$SRC/derived/clingen_gdv.tsv"   # check_sa_db.sh reads it
            FORCE_SAVE="$FORCE"; FORCE=1
            build omim omim oga "$SRC/derived/clingen_gdv.tsv" "${cg[0]}"
            FORCE="$FORCE_SAVE"
        fi
    fi

    RM="$SRC/rmsk.txt.gz"
    if up_to_date repeatmasker.osi "$RM" --name repeatmasker; then
        echo "[skip] repeatmasker.osi (current)" >&2; SKIPPED="$SKIPPED repeatmasker"
    else
        [ -s "$RM" ] || die "repeatmasker: missing $RM"
        python3 "$SA_SOURCES/repeatmasker_to_bed.py" "$RM" > "$WORK/repeatmasker.bed" \
            2> "$LOGDIR/repeatmasker_to_bed.log" || { cat "$LOGDIR/repeatmasker_to_bed.log" >&2; die "repeatmasker_to_bed.py failed"; }
        cat "$LOGDIR/repeatmasker_to_bed.log" >&2
        # --name repeatmasker is load-bearing: the classifier finds BP3's track by a key containing "repeat".
        FORCE_SAVE="$FORCE"; FORCE=1
        build custom_bed repeatmasker osi "$WORK/repeatmasker.bed" "$RM" --name repeatmasker
        FORCE="$FORCE_SAVE"
    fi
else
    c="$CHROM"
    G="$SRC/gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz"
    build gnomad "gnomad_chr${c}" osa2 "$G" "$G" --chunk-bits 16

    # SpliceAI and PhyloP are distilled from the same gnomAD file (its
    # spliceai_ds_max and phylop INFO fields). --url points the extractor at the
    # local .bgz; with no --regions it reads the whole chromosome via tabix.
    if up_to_date "spliceai_chr${c}.osa2" "$G" --chunk-bits 16 && \
       up_to_date "phylop_chr${c}.osa2" "$G" --chunk-bits 16; then
        echo "[skip] spliceai_chr${c}.osa2 phylop_chr${c}.osa2 (current)" >&2
        SKIPPED="$SKIPPED spliceai_chr${c} phylop_chr${c}"
    else
        [ -s "$G" ] && [ -s "$G.tbi" ] || die "chr$c: need $G and its .tbi for extract_gnomad_scores.py"
        echo "[extract] spliceai/phylop chr$c <- $G" >&2
        python3 "$SCRIPTS/extract_gnomad_scores.py" --chrom "$c" --url "$G" --out-dir "$WORK/extract" \
            2> "$LOGDIR/extract_gnomad_scores_chr${c}.log" || {
            cat "$LOGDIR/extract_gnomad_scores_chr${c}.log" >&2; die "extract_gnomad_scores.py chr$c failed"; }
        cat "$LOGDIR/extract_gnomad_scores_chr${c}.log" >&2
        FORCE_SAVE="$FORCE"; FORCE=1
        build spliceai "spliceai_chr${c}" osa2 "$WORK/extract/spliceai_chr${c}.vcf" "$G" --chunk-bits 16
        build phylop "phylop_chr${c}" osa2 "$WORK/extract/phylop_chr${c}.tsv" "$G" --chunk-bits 16
        FORCE="$FORCE_SAVE"
    fi

    if [ "$c" = Y ]; then
        echo "[skip] revel_chrY: REVEL v1.3 has no chrY records, and sa-build refuses an empty database" >&2
    else
        R="$SRC/revel/revel_chr${c}.csv"
        build revel "revel_chr${c}" osa2 "$R" "$R"
    fi

    if [ "$DELETE_GNOMAD" = 1 ]; then
        if is_done "gnomad_chr${c}.osa2" && is_done "spliceai_chr${c}.osa2" && is_done "phylop_chr${c}.osa2"; then
            [ -s "$G.sha256" ] || ident "$G" > /dev/null    # keep its identity for later skips
            echo "[clean] --delete-gnomad-src: removing $G(.tbi)" >&2
            rm -f "$G" "$G.tbi"
        else
            echo "[clean] NOT deleting $G: not every database built from it has a .done" >&2
        fi
    fi
fi

echo "[done] mode=$MODE${CHROM:+ chr$CHROM} built:${BUILT:- none} skipped:${SKIPPED:- none}" >&2
echo "[done] verify the whole stack with: acmg_fastvep/check_sa_db.sh" >&2
