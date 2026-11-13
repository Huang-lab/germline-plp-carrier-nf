#!/usr/bin/env bash
# Download every upstream file the fastVEP v0.5.0 ACMG reference stack is built
# from, into $FV_HOME/sa_src, and check each one before anything builds on it.
#
# Run on a Minerva LOGIN node (compute nodes have no internet). Resumable:
# curl -C - continues a partial download, and a file that is already complete
# and verified is not fetched again. Builds nothing — build_sa_dbs.sh does that.
#
# Why the checks: `fastvep sa-build` is a converter, not a downloader. Handed a
# truncated file it can succeed and write a database holding part of the
# source, and an ACMG run against that database classifies without the missing
# evidence and without a word (fastVEP docs/ACMG_SETUP.md, "The failure mode to
# design against"). So every file is checked here: md5 for ClinVar, gzip -t on
# the small gzips, tabix -l plus the BGZF end-of-file block for gnomAD
# (--full-check also runs bgzip -t, which reads all ~200 GB).
#
# Sources (the stack of fastVEP docs/ACMG_SETUP.md at v0.5.0):
#   ClinVar weekly VCF + .md5          clinvar_<YYYYMMDD>.vcf.gz       -> clinvar
#   ClinVar variant_summary.txt.gz     variant_summary.txt.gz          -> clinvar_protein
#   gnomAD v4.1 EXOMES chr1-22,X,Y     gnomad.exomes.v4.1.sites.chrN.vcf.bgz(.tbi)
#                                                                      -> gnomad, spliceai, phylop
#   gnomAD v4.1 constraint             gnomad.v4.1.constraint_metrics.tsv -> gnomad_genes
#   REVEL v1.3                         revel/revel_chrN.csv (split)    -> revel
#   ClinGen Gene-Disease Validity      clingen_gene_validity_<YYYYMMDD>.csv -> omim (json key)
#   UCSC RepeatMasker                  rmsk.txt.gz                     -> repeatmasker
#   (opt-in) OMIM genemap2.txt         genemap2.txt, only with --omim-genemap2 (licensed)
#
# Writes sa_src/SOURCES.tsv: file, url, date (server Last-Modified via curl -R,
# else download time), size, sha256, note. sha256 is cached per file in
# <file>.sha256 so a re-run does not re-hash 200 GB.
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: fetch_sa_sources.sh [options]      (needs FV_HOME, e.g. `source $FV_HOME/env.sh`)

Options:
  --clinvar-date <YYYYMMDD>  pin a dated weekly ClinVar release (recommended). Default:
                             the moving clinvar.vcf.gz, renamed after its ##fileDate.
  --chroms "<list>"          gnomAD chromosomes (default: "1 2 ... 22 X Y")
  --skip-gnomad              do not download gnomAD (e.g. already fetched)
  --full-check               also bgzip -t every gnomAD file (reads all of it)
  --omim-genemap2 <path>     copy a LICENSED OMIM genemap2.txt into sa_src (opt-in;
                             build_sa_dbs.sh uses it only with --use-omim)
  -h, --help
USAGE
    exit 2
}

: "${FV_HOME:?set FV_HOME (source \$FV_HOME/env.sh written by setup_fastvep.sh)}"
SRC="$FV_HOME/sa_src"
CLINVAR_DATE=""; CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
SKIP_GNOMAD=0; FULL_CHECK=0; OMIM=""
while [ $# -gt 0 ]; do
    case "$1" in
        --clinvar-date)  CLINVAR_DATE="$2"; shift 2 ;;
        --chroms)        CHROMS="$2"; shift 2 ;;
        --skip-gnomad)   SKIP_GNOMAD=1; shift ;;
        --full-check)    FULL_CHECK=1; shift ;;
        --omim-genemap2) OMIM="$2"; shift 2 ;;
        -h|--help)       usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
if [ -n "$CLINVAR_DATE" ]; then
    case "$CLINVAR_DATE" in
        [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]) ;;
        *) echo "ERROR: --clinvar-date must be YYYYMMDD (got $CLINVAR_DATE)" >&2; exit 2 ;;
    esac
fi
for t in curl gzip unzip tabix; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: $t not on PATH (module load htslib for tabix/bgzip)" >&2; exit 3; }
done
mkdir -p "$SRC"

CLINVAR_BASE="https://ftp.ncbi.nlm.nih.gov/pub/clinvar"
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1"
REVEL_URL="https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip"
CLINGEN_URL="https://search.clinicalgenome.org/kb/gene-validity/download"
RMSK_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"

REG="$(mktemp)"; trap 'rm -f "$REG"' EXIT
register() { printf '%s\t%s\t%s\n' "$1" "$2" "${3:-}" >> "$REG"; }   # file url note

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}
md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1
    else md5 -q "$1"; fi
}
fetch() {   # url dest
    echo "[fetch] $1" >&2
    curl -fL -R -C - --retry 3 --retry-delay 10 -o "$2" "$1"
}
fetch_small() {   # url dest — small sidecar files (.md5, .tbi): always fetched whole
    echo "[fetch] $1" >&2
    rm -f "$2"
    curl -fL -R --retry 3 --retry-delay 10 -o "$2" "$1"
}
die() { echo "ERROR: $*" >&2; exit 4; }

# --- ClinVar allele VCF ---
if [ -n "$CLINVAR_DATE" ]; then
    CV_URL="$CLINVAR_BASE/vcf_GRCh38/weekly/clinvar_${CLINVAR_DATE}.vcf.gz"
    CV="$SRC/clinvar_${CLINVAR_DATE}.vcf.gz"
    if [ ! -s "$CV.md5.ok" ]; then
        fetch_small "$CV_URL.md5" "$CV.md5"
        fetch "$CV_URL" "$CV"
    fi
else
    CV_URL="$CLINVAR_BASE/vcf_GRCh38/clinvar.vcf.gz"
    # The moving file: download, then name it after the release date it carries.
    TMPCV="$SRC/clinvar.latest.vcf.gz"
    rm -f "$TMPCV" "$TMPCV.md5"
    fetch_small "$CV_URL.md5" "$TMPCV.md5"
    fetch "$CV_URL" "$TMPCV"
    d="$(gzip -dc "$TMPCV" 2>/dev/null | head -50 | sed -n 's/^##fileDate=//p' | tr -d '-' | head -1 || true)"
    case "$d" in [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]) ;; *) die "no ##fileDate in $CV_URL header" ;; esac
    CLINVAR_DATE="$d"
    CV="$SRC/clinvar_${CLINVAR_DATE}.vcf.gz"
    mv -f "$TMPCV" "$CV"; mv -f "$TMPCV.md5" "$CV.md5"
    echo "[fetch] latest ClinVar is ##fileDate $CLINVAR_DATE -> $CV" >&2
fi
if [ ! -s "$CV.md5.ok" ]; then
    want="$(cut -d' ' -f1 "$CV.md5")"; got="$(md5_of "$CV")"
    [ "$want" = "$got" ] || die "$CV md5 $got != published $want (truncated download? delete it and re-run)"
    gzip -t "$CV" || die "$CV fails gzip -t"
    # The PM2 frequency backstop lives in AF_EXAC/AF_TGP/AF_ESP; a stripped VCF lacks it.
    gzip -dc "$CV" | head -300 | grep -q '^##INFO=<ID=AF_EXAC' || die "$CV has no AF_EXAC INFO header"
    echo "$got" > "$CV.md5.ok"
fi
# Keep exactly one ClinVar VCF in sa_src, so the build cannot pick up two.
for old in "$SRC"/clinvar_*.vcf.gz; do
    [ "$old" = "$CV" ] || [ ! -e "$old" ] || { echo "[fetch] NOTE: older ClinVar $old left in place; build_sa_dbs.sh uses CLINVAR_DATE=$CLINVAR_DATE" >&2; }
done
echo "$CLINVAR_DATE" > "$SRC/CLINVAR_DATE"
register "$CV" "$CV_URL" "ClinVar ##fileDate=$CLINVAR_DATE; md5 verified"

# --- ClinVar variant_summary (clinvar_protein) ---
VS="$SRC/variant_summary.txt.gz"; VS_URL="$CLINVAR_BASE/tab_delimited/variant_summary.txt.gz"
if [ ! -s "$VS.ok" ]; then
    rm -f "$VS"
    fetch "$VS_URL" "$VS"
    if curl -fsL -o "$VS.md5" "$VS_URL.md5"; then
        [ "$(cut -d' ' -f1 "$VS.md5")" = "$(md5_of "$VS")" ] || die "$VS md5 mismatch"
    else
        echo "[fetch] WARN: no $VS_URL.md5; relying on gzip -t" >&2; rm -f "$VS.md5"
    fi
    gzip -t "$VS" || die "$VS fails gzip -t"
    date -u +%Y%m%d > "$VS.ok"
fi
register "$VS" "$VS_URL" "moving file fetched $(cat "$VS.ok"); pair with ClinVar $CLINVAR_DATE (NCBI keeps no weekly copy)"

# --- gnomAD v4.1 exomes, per chromosome ---
bgzf_eof_ok() {   # the 28-byte empty BGZF block every complete bgzip file ends with
    [ "$(tail -c 28 "$1" | od -An -tx1 | tr -d ' \n')" = "1f8b08040000000000ff0600424302001b0003000000000000000000" ]
}
for c in $CHROMS; do
    f="$SRC/gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz"
    url="$GNOMAD_BASE/vcf/exomes/gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz"
    if [ "$SKIP_GNOMAD" = 1 ]; then
        if [ -s "$f" ] || [ -s "$f.sha256" ]; then
            register "$f" "$url" "gnomAD v4.1 exomes (verified earlier: $( [ -e "$f.ok" ] && echo yes || echo NO))"
        fi
        continue
    fi
    if [ ! -e "$f.ok" ]; then
        fetch_small "$url.tbi" "$f.tbi"
        fetch "$url" "$f"
        bgzf_eof_ok "$f" || die "$f does not end in a BGZF EOF block (truncated)"
        tabix -l "$f" | grep -qx "chr${c}" || die "tabix -l $f does not list chr${c}"
        if [ "$FULL_CHECK" = 1 ]; then
            command -v bgzip >/dev/null 2>&1 || die "--full-check needs bgzip"
            bgzip -t "$f" || die "$f fails bgzip -t"
        fi
        date -u +%Y%m%d > "$f.ok"
    fi
    register "$f" "$url" "gnomAD v4.1 exomes; EOF+tabix verified$( [ "$FULL_CHECK" = 1 ] && echo '+bgzip -t')"
    if [ -e "$f.tbi" ] || [ -s "$f.tbi.sha256" ]; then register "$f.tbi" "$url.tbi" ""; fi
done

# --- gnomAD v4.1 gene constraint ---
GC="$SRC/gnomad.v4.1.constraint_metrics.tsv"; GC_URL="$GNOMAD_BASE/constraint/gnomad.v4.1.constraint_metrics.tsv"
if [ ! -e "$GC.ok" ]; then
    rm -f "$GC"; fetch "$GC_URL" "$GC"
    n="$(wc -l < "$GC" | tr -d ' ')"
    [ "$n" -gt 10000 ] || die "$GC has only $n lines"
    head -1 "$GC" | grep -q 'gene' || die "$GC header has no gene column"
    date -u +%Y%m%d > "$GC.ok"
fi
register "$GC" "$GC_URL" "gnomAD v4.1 constraint"

# --- REVEL v1.3: zip -> one CSV per chromosome (the build reads revel/revel_chrN.csv) ---
RZ="$SRC/revel-v1.3_all_chromosomes.zip"
if [ ! -e "$SRC/revel/.ok" ]; then
    fetch "$REVEL_URL" "$RZ"
    unzip -tq "$RZ" >/dev/null || die "$RZ fails unzip -t"
    rm -rf "$SRC/revel.partial"; mkdir -p "$SRC/revel.partial"
    # chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
    unzip -p "$RZ" revel_with_transcript_ids | awk -F, -v d="$SRC/revel.partial" '
        NR == 1 { hdr = $0; next }
        { f = d "/revel_chr" $1 ".csv"
          if (!(f in seen)) { print hdr > f; seen[f] = 1 }
          print >> f; n[$1]++ }
        END { for (c in n) printf "%s\t%d\n", c, n[c] > (d "/COUNTS.tsv") }'
    [ -s "$SRC/revel.partial/revel_chr1.csv" ] || die "REVEL split produced no chr1"
    rm -rf "$SRC/revel"; mv "$SRC/revel.partial" "$SRC/revel"
    date -u +%Y%m%d > "$SRC/revel/.ok"
fi
register "$RZ" "$REVEL_URL" "REVEL v1.3; split to revel/revel_chrN.csv (counts: revel/COUNTS.tsv)"

# --- ClinGen Gene-Disease Validity (the `omim` json key's default content) ---
shopt -s nullglob; cg=("$SRC"/clingen_gene_validity_*.csv); shopt -u nullglob
if [ ${#cg[@]} -eq 0 ]; then
    CG="$SRC/clingen_gene_validity_$(date -u +%Y%m%d).csv"
    rm -f "$CG.partial"
    echo "[fetch] $CLINGEN_URL" >&2
    curl -fL --retry 3 --retry-delay 10 -o "$CG.partial" "$CLINGEN_URL"
    grep -q 'GENE SYMBOL' "$CG.partial" || die "ClinGen download has no 'GENE SYMBOL' header row (got an error page?)"
    n="$(wc -l < "$CG.partial" | tr -d ' ')"
    [ "$n" -gt 1000 ] || die "ClinGen download has only $n lines"
    mv -f "$CG.partial" "$CG"
elif [ ${#cg[@]} -gt 1 ]; then
    die "more than one clingen_gene_validity_*.csv in $SRC; keep the one to build from"
else
    CG="${cg[0]}"
fi
register "$CG" "$CLINGEN_URL" "ClinGen GDV export, dated by download day"

# --- UCSC RepeatMasker ---
RM="$SRC/rmsk.txt.gz"
if [ ! -e "$RM.ok" ]; then
    rm -f "$RM"; fetch "$RMSK_URL" "$RM"
    gzip -t "$RM" || die "$RM fails gzip -t"
    date -u +%Y%m%d > "$RM.ok"
fi
register "$RM" "$RMSK_URL" "UCSC hg38 rmsk"

# --- OMIM genemap2 (opt-in, licensed; never downloaded) ---
if [ -n "$OMIM" ]; then
    [ -s "$OMIM" ] || die "--omim-genemap2: $OMIM missing or empty"
    cp -f "$OMIM" "$SRC/genemap2.txt"
    register "$SRC/genemap2.txt" "file://$OMIM" "OMIM genemap2 (licensed; used only with build_sa_dbs.sh --use-omim)"
fi

# --- SOURCES.tsv ---
OUT="$SRC/SOURCES.tsv"
{
    printf 'file\turl\tdate\tsize\tsha256\tnote\n'
    while IFS=$'\t' read -r f url note; do
        if [ ! -e "$f" ] && [ -s "$f.sha256" ]; then
            # gnomAD source removed by build_sa_dbs.sh --delete-gnomad-src after its build
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${f#"$SRC"/}" "$url" "-" "deleted" "$(cat "$f.sha256")" \
                "$note; deleted after build"
            continue
        fi
        if [ ! -s "$f.sha256" ] || [ "$f.sha256" -ot "$f" ]; then
            echo "[fetch] sha256 $f" >&2
            sha256_of "$f" > "$f.sha256"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${f#"$SRC"/}" "$url" "$(date -r "$f" +%Y-%m-%d)" \
            "$(wc -c < "$f" | tr -d ' ')" "$(cat "$f.sha256")" "$note"
    done < "$REG"
} > "$OUT.partial"
mv -f "$OUT.partial" "$OUT"
echo "[fetch] done. $(($(wc -l < "$OUT") - 1)) files recorded in $OUT" >&2
