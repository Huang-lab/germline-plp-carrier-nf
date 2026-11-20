#!/usr/bin/env bash
# Standalone fastVEP runner for a SINGLE QC'd VCF (Huang-lab/fastVEP v0.5.0).
#
# NOT part of the Nextflow pipeline. Runs after the pipeline's QC step, on the
# QC'd VCFs it publishes (results/<run>/norm_qc/<chunk>.norm.vcf.gz), so fastVEP
# annotates exactly the same variants as the ClinVar path — without touching the
# working Nextflow pipeline.
#
# Two modes:
#   * plain annotation (default)  — transcript cache + FASTA, no supplementary
#     databases. Emits <name>.fastvep.vcf.gz.
#   * ACMG (--acmg)               — additionally reads the sa_db built by
#     build_sa_dbs.sh and verified by check_sa_db.sh, and runs with
#     --flag-pick-allele-gene: every transcript's CSQ block is kept and fastVEP
#     marks its own choice per (allele, gene) with PICK=1, which the carriers
#     step (build_acmg_carriers_pick.py, fastvep_pick mode) reads instead of
#     re-deriving a pick. Emits <name>.fastvep.vcf.gz AND <name>.acmg_plp.tsv.
#
# Everything comes from env.sh (setup_fastvep.sh) — no PATH lookup, no
# "first file in the directory" guessing:
#   FV_ENV=<$FV_HOME/env.sh>   sourced if set; or export these yourself:
#   FV_BIN FV_BIN_SHA256 FV_GFF3 FV_FASTA FV_CACHE [FV_SA for --acmg]
# Command-line options override the environment.
#
# A chunk FAILS (nothing published) when fastVEP's stderr says the reference
# stack was not what we think it is — fastVEP itself still exits 0 then:
#   'incomplete reference stack'           a source of the nine is not loaded
#   'was not loaded'                       a database file is empty/unopenable
#   'hold nothing of chromosomes'          gnomAD lacks a chromosome of the gene model
#   'answered by more than one database'   two databases answer for one source key
#   'could not be read for'                an SA lookup hit an I/O error: those variants
#                                          carry no annotation (PM2 then fires falsely)
#   'could not be read and were not annotated'  VCF records dropped from the output
#   'not a contig of the gene models' / 'are on contigs that the gene models and the
#   reference do not have'                 records silently called intergenic
#   'holds databases built for different assemblies'  sa_db mixes assemblies
#                                          (check_sa_assemblies, printed only on a conflict)
#   'have a REF allele that is not the reference'      REF disagrees with --fasta
#                                          (printed only when >= 1 record with plain-base
#                                          REF differs; N/ambiguity codes are not counted;
#                                          our QC'd VCFs are bcftools-normalised against the
#                                          same GRCh38, so any count means another assembly/FASTA)
# (strings from fastvep-annotate/src/lib.rs and fastvep-cli/src/pipeline/annotate.rs
# at v0.5.0). fastVEP's stderr is kept as <name>.fastvep.stderr.log either way.
#
# One writer per chunk: a lock on <outdir>/.<name>.lock (fv_lock.sh), and
# partial files named per job, so a second job on the same chunk is refused
# instead of publishing the first job's half-written file.
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: run_fastvep.sh -i <qc.vcf[.gz]> -o <outdir> [options]

Required:
  -i, --input <vcf>       QC'd VCF (.vcf or .vcf.gz), Ensembl contigs
  -o, --outdir <dir>      output directory (created if absent)
  plus FV_BIN FV_BIN_SHA256 FV_GFF3 FV_FASTA FV_CACHE (env, or FV_ENV=<env.sh>)

Options (each overrides its env variable):
  --fastvep <path>        fastvep binary                       ($FV_BIN)
  --gff3 <gff3>           GFF3 the transcript cache was built from, recorded only ($FV_GFF3)
  --fasta <fa>            reference FASTA                      ($FV_FASTA)
  --transcript-cache <f>  fastVEP transcript cache, always used ($FV_CACHE)
  --sa-dir <dir>          supplementary-annotation DB dir (--acmg) ($FV_SA)
  --acmg                  ACMG-AMP classification, with --flag-pick-allele-gene
  --acmg-config <toml>    ACMG threshold overrides
  --hgvs                  emit HGVSc/HGVSp (--acmg turns it on anyway)
  --pick                  plain mode only: one consequence per variant
  --threads <N>           worker threads (exported as RAYON_NUM_THREADS)
  --expect-provenance     print the provenance key fields this run WOULD record and
                          exit without running (run_fastvep_batch.sh's skip check)
  -h, --help

Outputs:
  <outdir>/<name>.fastvep.vcf.gz        annotated VCF
  <outdir>/<name>.acmg_plp.tsv          per-variant ACMG table (--acmg only)
  <outdir>/<name>.fastvep.provenance    binary/SA hashes this chunk was made with
  <outdir>/<name>.fastvep.stderr.log    fastVEP's stderr
USAGE
    exit 2
}

if [ -n "${FV_ENV:-}" ]; then
    # shellcheck disable=SC1090
    . "$FV_ENV"
fi
FASTVEP="${FV_BIN:-}"; GFF3="${FV_GFF3:-}"; FASTA="${FV_FASTA:-}"
TRANSCRIPT_CACHE="${FV_CACHE:-}"; SA_DIR="${FV_SA:-}"
IN=""; OUTDIR=""; ACMG_CONFIG=""; THREADS=""
DO_ACMG=0; DO_HGVS=0; DO_PICK=0; EXPECT_ONLY=0
while [ $# -gt 0 ]; do
    case "$1" in
        -i|--input)          IN="$2"; shift 2 ;;
        -o|--outdir)         OUTDIR="$2"; shift 2 ;;
        --gff3)              GFF3="$2"; shift 2 ;;
        --fasta)             FASTA="$2"; shift 2 ;;
        --sa-dir)            SA_DIR="$2"; shift 2 ;;
        --acmg-config)       ACMG_CONFIG="$2"; shift 2 ;;
        --transcript-cache)  TRANSCRIPT_CACHE="$2"; shift 2 ;;
        --threads)           THREADS="$2"; shift 2 ;;
        --fastvep)           FASTVEP="$2"; shift 2 ;;
        --acmg)              DO_ACMG=1; shift ;;
        --hgvs)              DO_HGVS=1; shift ;;
        --pick)              DO_PICK=1; shift ;;
        --expect-provenance) EXPECT_ONLY=1; shift ;;
        -h|--help)           usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}
file_stat() { printf '%s|%s|%s' "$1" "$(wc -c < "$1" | tr -d ' ')" "$(date -r "$1" +%s)"; }

# --- validation ---
[ -n "$IN" ] && [ -n "$OUTDIR" ] || usage
missing=""
for v in FASTVEP:FV_BIN GFF3:FV_GFF3 FASTA:FV_FASTA TRANSCRIPT_CACHE:FV_CACHE; do
    eval "val=\${${v%%:*}}"
    [ -n "$val" ] || missing="$missing ${v#*:}"
done
[ -n "${FV_BIN_SHA256:-}" ] || missing="$missing FV_BIN_SHA256"
[ -z "$missing" ] || { echo "ERROR: not set:$missing — source \$FV_HOME/env.sh (or set FV_ENV)" >&2; exit 3; }
[ -x "$FASTVEP" ] || { echo "ERROR: fastvep binary not executable: $FASTVEP" >&2; exit 3; }
for f in "$IN" "$GFF3" "$FASTA" "$TRANSCRIPT_CACHE"; do [ -s "$f" ] || { echo "ERROR: missing/empty: $f" >&2; exit 3; }; done
BIN_SHA="$(sha256_of "$FASTVEP")"
[ "$BIN_SHA" = "$FV_BIN_SHA256" ] || {
    echo "ERROR: sha256($FASTVEP)=$BIN_SHA != FV_BIN_SHA256=$FV_BIN_SHA256 — not the pinned binary" >&2; exit 3; }
SA_SHA="none"
if [ "$DO_ACMG" = 1 ]; then
    [ -n "$SA_DIR" ] || { echo "ERROR: --acmg requires --sa-dir / FV_SA" >&2; exit 3; }
    [ -d "$SA_DIR" ] || { echo "ERROR: --sa-dir not a directory: $SA_DIR" >&2; exit 3; }
    [ -s "$SA_DIR/SA_DB_SHA256.txt" ] || {
        echo "ERROR: $SA_DIR has no SA_DB_SHA256.txt — run check_sa_db.sh on it first" >&2; exit 3; }
    SA_SHA="$(sha256_of "$SA_DIR/SA_DB_SHA256.txt")"
    # --flag-pick-allele-gene and --pick are one clap group in fastVEP: exclusive.
    [ "$DO_PICK" = 0 ] || { echo "ERROR: --pick cannot be combined with --acmg (which flags picks instead)" >&2; exit 3; }
fi
[ -n "$THREADS" ] && export RAYON_NUM_THREADS="$THREADS"

HERE="$(cd "$(dirname "$0")" && pwd)"
name="$(basename "$IN")"; name="${name%.gz}"; name="${name%.vcf}"; name="${name%.norm}"

# --- the fastVEP flags (everything but input/output), and the provenance key ---
# The transcript cache is the transcript source; fastVEP ignores --gff3 when a
# cache loads (and warns), so --gff3 is not passed. FV_GFF3 is recorded only.
flags=(--output-format vcf --transcript-cache "$TRANSCRIPT_CACHE" --fasta "$FASTA" --no-progress)
[ "$DO_HGVS" = 1 ]           && flags+=(--hgvs)
[ "$DO_PICK" = 1 ]           && flags+=(--pick)
if [ "$DO_ACMG" = 1 ]; then
    flags+=(--sa-dir "$SA_DIR" --acmg --flag-pick-allele-gene)
    [ -n "$ACMG_CONFIG" ] && flags+=(--acmg-config "$ACMG_CONFIG")
fi
ACMG_CONFIG_SHA="none"
if [ -n "$ACMG_CONFIG" ]; then
    [ -s "$ACMG_CONFIG" ] || { echo "ERROR: missing/empty --acmg-config: $ACMG_CONFIG" >&2; exit 3; }
    ACMG_CONFIG_SHA="$(sha256_of "$ACMG_CONFIG")"
fi
# The cache is several GB: run_fastvep_batch.sh hashes it once and passes the
# digest with the path|size|mtime it was taken at; anything else is re-hashed.
if [ -n "${FV_CACHE_SHA256:-}" ] && [ "${FV_CACHE_STAT:-}" = "$(file_stat "$TRANSCRIPT_CACHE")" ]; then
    CACHE_SHA="$FV_CACHE_SHA256"
else
    CACHE_SHA="$(sha256_of "$TRANSCRIPT_CACHE")"
fi
# Every field here must match for run_fastvep_batch.sh to skip a finished chunk.
provenance_key() {
    echo "mode             $([ "$DO_ACMG" = 1 ] && echo acmg || echo annotate)"
    echo "fastvep_sha256   $BIN_SHA"
    echo "sa_db_sha256     $SA_SHA"
    echo "input_sha256     $(sha256_of "$IN")"
    echo "transcript_cache_sha256 $CACHE_SHA"
    echo "fasta            $FASTA"
    echo "acmg_config_sha256 $ACMG_CONFIG_SHA"
    echo "fastvep_flags    ${flags[*]}"
}
if [ "$EXPECT_ONLY" = 1 ]; then
    provenance_key
    exit 0
fi

mkdir -p "$OUTDIR"
# Write the annotated VCF gzip-compressed. fastVEP streams to stdout (--output -)
# and we gzip on the fly, so the large uncompressed VCF NEVER lands on disk —
# keeps peak disk low and makes writing into space-tight result dirs safe.
# The parser reads .gz transparently.
VEPOUT="$OUTDIR/${name}.fastvep.vcf.gz"
PROVOUT="$OUTDIR/${name}.fastvep.provenance"
ERRLOG="$OUTDIR/${name}.fastvep.stderr.log"

# Outputs are published by rename, never streamed into their final path.
#
# The redirect `... | gzip -c > "$VEPOUT"` creates $VEPOUT the instant the run
# starts, so a chunk killed at the LSF walltime leaves a non-empty file holding
# a prefix of the records. Worse, that prefix is a *valid* gzip stream — `gzip
# -t` passes on it, because gzip closed cleanly when its stdin went away — so
# nothing downstream can tell it from a complete annotation. run_fastvep_batch.sh
# then reads `[ -s "$VEPOUT" ]` as "already done" and skips the chunk on the
# resubmit that was supposed to finish it. That is the batch3 dropped-chunk
# incident with the manifest reporting success.
#
# Streaming to a temp file and renaming on completion makes the final path mean
# what the batch driver assumes it means: present iff complete. Same lesson as
# Huang-lab/fastVEP#88, one layer up. The temp names carry the job id and pid,
# and the chunk lock below makes this job the only writer.
# shellcheck source=fv_lock.sh
. "$HERE/fv_lock.sh"
fv_lock_acquire "$OUTDIR/.${name}.lock" 9 "chunk $name in $OUTDIR" || exit 7
TAG="${LSB_JOBID:-local}.$$"
VEPTMP="$OUTDIR/.${name}.fastvep.vcf.gz.partial.$TAG"
ACMGTMP="$OUTDIR/.${name}.acmg_plp.tsv.partial.$TAG"
PROVTMP="$OUTDIR/.${name}.fastvep.provenance.partial.$TAG"
ERRTMP="$OUTDIR/.${name}.fastvep.stderr.log.partial.$TAG"
trap 'rm -f "$VEPTMP" "$ACMGTMP" "$PROVTMP" "$ERRTMP"; fv_lock_release' EXIT
trap 'exit 143' TERM INT
# Holding the lock, any other partial of this chunk is from a dead job: clear it.
# The provenance goes first on a re-run: outputs without it are never skipped.
rm -f "$OUTDIR/.${name}".*.partial* "$PROVOUT"

cmd=("$FASTVEP" annotate --input "$IN" --output - "${flags[@]}")
echo "[fastvep] ${cmd[*]} 2> $ERRLOG | gzip > $VEPOUT" >&2
if ! "${cmd[@]}" 2> "$ERRTMP" | gzip -c > "$VEPTMP"; then
    mv -f "$ERRTMP" "$ERRLOG"
    echo "ERROR: fastvep failed for $name (stderr: $ERRLOG)" >&2
    exit 1
fi
mv -f "$ERRTMP" "$ERRLOG"

# --- refuse a run whose reference stack was not the one we think ---
GATED='incomplete reference stack|was not loaded|hold nothing of chromosomes|answered by more than one database'
GATED="$GATED|could not be read for|could not be read and were not annotated"
GATED="$GATED|not a contig of the gene models|are on contigs that the gene models and the reference do not have"
GATED="$GATED|holds databases built for different assemblies|have a REF allele that is not the reference"
if grep -E "$GATED" "$ERRLOG" >&2; then
    echo "ERROR: fastVEP reported a reference-stack problem (lines above, full log $ERRLOG);" \
         "not publishing $name" >&2
    exit 5
fi

# --- ACMG post-parse (only in --acmg mode) ---
#
# Both artifacts are produced before either is published, so the pair the batch
# driver tests for is never half-present: a parse that dies does not leave a
# complete VCF next to a truncated TSV.
if [ "$DO_ACMG" = 1 ]; then
    echo "[parse] (pending) -> $OUTDIR/${name}.acmg_plp.tsv" >&2
    python3 "$HERE/parse_fastvep_acmg.py" --vcf "$VEPTMP" --out "$ACMGTMP"
fi

{
    provenance_key
    echo "chunk            $name"
    echo "input            $IN"
    echo "fastvep_version  $("$FASTVEP" --version 2>/dev/null || echo unknown)"
    echo "sa_dir           ${SA_DIR:-none}"
    echo "transcript_cache $TRANSCRIPT_CACHE"
    echo "command          ${cmd[*]}"
    echo "date_utc         $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "lsb_jobid        ${LSB_JOBID:-}"
} > "$PROVTMP"

if [ "$DO_ACMG" = 1 ]; then
    mv -f "$VEPTMP" "$VEPOUT"
    mv -f "$ACMGTMP" "$OUTDIR/${name}.acmg_plp.tsv"
    mv -f "$PROVTMP" "$PROVOUT"
    echo "[done] $VEPOUT  $OUTDIR/${name}.acmg_plp.tsv" >&2
else
    mv -f "$VEPTMP" "$VEPOUT"
    mv -f "$PROVTMP" "$PROVOUT"
    echo "[done] $VEPOUT" >&2
fi
