#!/usr/bin/env bash
# Batch driver: run fastVEP over EVERY QC'd VCF in a norm_qc/ directory.
#
# For each <chunk>.norm.vcf.gz it calls run_fastvep.sh (plain or --acmg). Runs
# locally in sequence by default, or submits one LSF job per chunk with --lsf.
# Writes a chunk->status manifest so a failed chunk is never silent (the
# lesson from the batch3 dropped-chunk incident). In --acmg mode, optionally
# concatenates the per-chunk *.acmg_plp.tsv into one table.
#
# Before any chunk runs: sha256(FV_BIN) must equal FV_BIN_SHA256 (env.sh), and
# in --acmg mode `check_sa_db.sh --quick` must pass on the sa_db. The batch
# writes <outdir>/fastvep_run_provenance.txt (binary, cache, FASTA, sa_db hash,
# pipeline git rev, args) and only skips a chunk whose own provenance matches
# the current binary + sa_db.
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: run_fastvep_batch.sh --in-dir <norm_qc/> -o <outdir> [options]

Required:
  --in-dir <dir>          directory of QC'd VCFs (*.norm.vcf.gz / *.vcf.gz)
  -o, --outdir <dir>      output directory
  plus FV_BIN FV_BIN_SHA256 FV_GFF3 FV_FASTA FV_CACHE [FV_SA] (env, or FV_ENV=<env.sh>)

Passthrough to run_fastvep.sh (each overrides its env variable):
  --gff3 <g>  --fasta <fa>  --transcript-cache <f>  --sa-dir <dir>  --fastvep <path>
  --hgvs  --pick  --acmg  --acmg-config <toml>  --threads <N>

Batch options:
  --lsf                   submit one bsub per chunk (else run locally, sequential)
  -P, --project <alloc>   LSF allocation (default: $MINERVA_ALLOCATION)
  --queue <q>             LSF queue (default: premium)
  --walltime <hh:mm>      LSF walltime (default: 4:00)
  --concat                after (local) runs, concatenate *.acmg_plp.tsv into one
  --force                 re-run chunks whose existing output was made with a different
                          binary / sa_db (or has no provenance) instead of refusing
  -h, --help
USAGE
    exit 2
}

if [ -n "${FV_ENV:-}" ]; then
    # shellcheck disable=SC1090
    . "$FV_ENV"
fi
HERE="$(cd "$(dirname "$0")" && pwd)"
INDIR=""; OUTDIR=""; ACMG_CONFIG=""; THREADS=""
FASTVEP="${FV_BIN:-}"; GFF3="${FV_GFF3:-}"; FASTA="${FV_FASTA:-}"
TRANSCRIPT_CACHE="${FV_CACHE:-}"; SA_DIR="${FV_SA:-}"
DO_ACMG=0; DO_HGVS=0; DO_PICK=0; USE_LSF=0; DO_CONCAT=0; FORCE=0
PROJECT="${MINERVA_ALLOCATION:-}"; QUEUE="premium"; WALL="4:00"
ORIG_ARGS="$*"
while [ $# -gt 0 ]; do
    case "$1" in
        --in-dir)             INDIR="$2"; shift 2 ;;
        -o|--outdir)          OUTDIR="$2"; shift 2 ;;
        --gff3)               GFF3="$2"; shift 2 ;;
        --fasta)              FASTA="$2"; shift 2 ;;
        --sa-dir)             SA_DIR="$2"; shift 2 ;;
        --acmg-config)        ACMG_CONFIG="$2"; shift 2 ;;
        --transcript-cache)   TRANSCRIPT_CACHE="$2"; shift 2 ;;
        --threads)            THREADS="$2"; shift 2 ;;
        --fastvep)            FASTVEP="$2"; shift 2 ;;
        --acmg)               DO_ACMG=1; shift ;;
        --hgvs)               DO_HGVS=1; shift ;;
        --pick)               DO_PICK=1; shift ;;
        --lsf)                USE_LSF=1; shift ;;
        -P|--project)         PROJECT="$2"; shift 2 ;;
        --queue)              QUEUE="$2"; shift 2 ;;
        --walltime)           WALL="$2"; shift 2 ;;
        --concat)             DO_CONCAT=1; shift ;;
        --force)              FORCE=1; shift ;;
        -h|--help)            usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$INDIR" ] && [ -n "$OUTDIR" ] || usage
[ -d "$INDIR" ] || { echo "ERROR: --in-dir not a directory: $INDIR" >&2; exit 3; }
if [ "$USE_LSF" = 1 ] && [ -z "$PROJECT" ]; then
    echo "ERROR: --lsf needs -P/--project or \$MINERVA_ALLOCATION" >&2; exit 3
fi
for v in FASTVEP:FV_BIN GFF3:FV_GFF3 FASTA:FV_FASTA TRANSCRIPT_CACHE:FV_CACHE; do
    eval "val=\${${v%%:*}}"
    [ -n "$val" ] || { echo "ERROR: ${v#*:} not set — source \$FV_HOME/env.sh (or set FV_ENV)" >&2; exit 3; }
done
: "${FV_BIN_SHA256:?FV_BIN_SHA256 not set — source \$FV_HOME/env.sh (or set FV_ENV)}"
export FV_BIN_SHA256

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

# --- pre-flight: the pinned binary, and an sa_db that still matches its manifest ---
[ -x "$FASTVEP" ] || { echo "ERROR: fastvep binary not executable: $FASTVEP" >&2; exit 3; }
BIN_SHA="$(sha256_of "$FASTVEP")"
if [ "$BIN_SHA" != "$FV_BIN_SHA256" ]; then
    echo "ERROR: sha256($FASTVEP)=$BIN_SHA != FV_BIN_SHA256=$FV_BIN_SHA256 — not the pinned binary" >&2
    exit 3
fi
SA_SHA="none"
if [ "$DO_ACMG" = 1 ]; then
    [ -n "$SA_DIR" ] || { echo "ERROR: --acmg requires --sa-dir / FV_SA" >&2; exit 3; }
    "$HERE/check_sa_db.sh" --quick --sa-dir "$SA_DIR" > /dev/null || {
        echo "ERROR: check_sa_db.sh --quick failed for $SA_DIR; not running the batch" >&2; exit 3; }
    SA_SHA="$(sha256_of "$SA_DIR/SA_DB_SHA256.txt")"
fi
MODE_NAME="$([ "$DO_ACMG" = 1 ] && echo acmg || echo annotate)"
mkdir -p "$OUTDIR"

# One batch driver per outdir at a time (the manifest and run provenance are
# shared files); each chunk is additionally locked by run_fastvep.sh itself.
# shellcheck source=fv_lock.sh
. "$HERE/fv_lock.sh"
fv_lock_acquire "$OUTDIR/.batch.lock" 8 "batch output dir $OUTDIR" || exit 7
trap 'fv_lock_release' EXIT
trap 'exit 143' TERM INT

# The transcript cache is hashed once here and handed to every chunk with the
# path|size|mtime it was taken at (run_fastvep.sh re-hashes if that differs).
FV_CACHE_SHA256="$(sha256_of "$TRANSCRIPT_CACHE")"
FV_CACHE_STAT="$(printf '%s|%s|%s' "$TRANSCRIPT_CACHE" "$(wc -c < "$TRANSCRIPT_CACHE" | tr -d ' ')" "$(date -r "$TRANSCRIPT_CACHE" +%s)")"
export FV_CACHE_SHA256 FV_CACHE_STAT

# common run_fastvep.sh args shared by every chunk
common=(--fastvep "$FASTVEP" --gff3 "$GFF3" --fasta "$FASTA" --transcript-cache "$TRANSCRIPT_CACHE")
[ "$DO_HGVS" = 1 ]         && common+=(--hgvs)
[ "$DO_PICK" = 1 ]         && common+=(--pick)
[ "$DO_ACMG" = 1 ]         && common+=(--acmg --sa-dir "$SA_DIR")
[ -n "$ACMG_CONFIG" ]      && common+=(--acmg-config "$ACMG_CONFIG")
[ -n "$THREADS" ]          && common+=(--threads "$THREADS")

# Run-level provenance: what every chunk of this batch is annotated with.
RUNPROV="$OUTDIR/fastvep_run_provenance.txt"
pipeline_rev="$(git -C "$HERE" rev-parse HEAD 2>/dev/null || echo unknown)"
pipeline_dirty="$( [ -n "$(git -C "$HERE" status --porcelain 2>/dev/null)" ] && echo yes || echo no)"
{
    echo "date_utc               $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "host                   $(hostname 2>/dev/null || echo unknown)"
    echo "lsb_jobid              ${LSB_JOBID:-}"
    echo "mode                   $MODE_NAME"
    echo "fastvep_bin            $FASTVEP"
    echo "fastvep_sha256         $BIN_SHA"
    echo "fastvep_version        $("$FASTVEP" --version 2>/dev/null || echo unknown)"
    echo "gff3                   $(basename "$GFF3")  ($GFF3)"
    echo "transcript_cache       $(basename "$TRANSCRIPT_CACHE")  ($TRANSCRIPT_CACHE)"
    echo "transcript_cache_sha256 $FV_CACHE_SHA256"
    echo "fasta                  $(basename "$FASTA")  ($FASTA)"
    echo "sa_dir                 ${SA_DIR:-none}"
    echo "sa_db_sha256           $SA_SHA   (sha256 of SA_DB_SHA256.txt)"
    echo "pipeline_rev           $pipeline_rev  (uncommitted changes: $pipeline_dirty)"
    echo "in_dir                 $INDIR"
    echo "args                   $ORIG_ARGS"
} > "$RUNPROV.partial"
mv -f "$RUNPROV.partial" "$RUNPROV"
echo "[batch] provenance: $RUNPROV" >&2

# collect inputs (norm.vcf.gz preferred; fall back to any vcf.gz / vcf)
shopt -s nullglob
inputs=("$INDIR"/*.norm.vcf.gz)
[ ${#inputs[@]} -eq 0 ] && inputs=("$INDIR"/*.vcf.gz)
[ ${#inputs[@]} -eq 0 ] && inputs=("$INDIR"/*.vcf)
shopt -u nullglob
[ ${#inputs[@]} -gt 0 ] || { echo "ERROR: no VCFs found in $INDIR" >&2; exit 3; }
echo "[batch] ${#inputs[@]} input VCF(s) from $INDIR" >&2

chunk_name() { local n; n="$(basename "$1")"; n="${n%.gz}"; n="${n%.vcf}"; echo "${n%.norm}"; }
# Value of a "field<spaces>value" line; awk reads to the end (no SIGPIPE).
kv() { awk -v k="$1" '$1 == k && !done { sub(/^[^ ]+ +/, ""); print; done = 1 }'; }
# The fields run_fastvep.sh --expect-provenance prints; all must match to skip.
KEY_FIELDS="mode fastvep_sha256 sa_db_sha256 input_sha256 transcript_cache_sha256 fasta acmg_config_sha256 fastvep_flags"

# Decide every chunk before running any. A chunk whose expected output(s)
# already exist is skipped only when its .fastvep.provenance says it was made
# from this input (sha256 of the chunk VCF), by this binary, against this
# sa_db and transcript cache (sha256), FASTA (path), acmg-config (sha256 or
# none), in this mode, with these fastVEP flags; output from anything else
# is refused (or redone with --force), so a re-run can never quietly mix
# annotations from two fastVEP versions or two sa_db builds in one batch.
#
# Presence still means "complete" only because run_fastvep.sh publishes by
# rename: before that, a chunk killed at the walltime left a truncated-but-valid
# .vcf.gz here and the resubmit meant to finish it skipped it instead —
# silently dropping every carrier in the missing tail.
todo=(); skip=(); refused=""
for vcf in "${inputs[@]}"; do
    name="$(chunk_name "$vcf")"
    vep_out="$OUTDIR/${name}.fastvep.vcf.gz"
    acmg_out="$OUTDIR/${name}.acmg_plp.tsv"
    prov="$OUTDIR/${name}.fastvep.provenance"
    for stale in "$OUTDIR/.${name}".*.partial*; do
        [ -e "$stale" ] && echo "[batch] note: stale partial from an interrupted run (cleared when the chunk reruns): $stale" >&2
    done
    have_out=0
    if [ -s "$vep_out" ]; then
        if [ "$DO_ACMG" = 1 ]; then [ -s "$acmg_out" ] && have_out=1; else have_out=1; fi
    fi
    if [ "$have_out" = 0 ]; then todo+=("$vcf"); continue; fi
    why=""
    if [ ! -s "$prov" ]; then
        why="no provenance file"
    else
        expect="$("$HERE/run_fastvep.sh" --expect-provenance -i "$vcf" -o "$OUTDIR" "${common[@]}")"
        for k in $KEY_FIELDS; do
            have="$(kv "$k" < "$prov")"; want="$(printf '%s\n' "$expect" | kv "$k")"
            if [ "$have" != "$want" ]; then
                why="$why${why:+; }$k differs (output: ${have:-<none>}; now: $want)"
            fi
        done
    fi
    if [ -z "$why" ]; then
        skip+=("$name")
    elif [ "$FORCE" = 1 ]; then
        echo "[batch] --force: redoing $name ($why)" >&2
        todo+=("$vcf")
    else
        refused="$refused
  $name: $why"
    fi
done
if [ -n "$refused" ]; then
    echo "[batch] ERROR: existing output in $OUTDIR was not made with this binary/sa_db/mode:$refused" >&2
    echo "[batch] Refusing to mix it into this batch. Use a fresh --outdir, or --force to redo those chunks." >&2
    exit 6
fi

MANIFEST="$OUTDIR/batch_manifest.tsv"
: > "$MANIFEST"
n_submitted=0; n_ok=0; n_fail=0; n_skipped=${#skip[@]}
for name in ${skip[@]+"${skip[@]}"}; do
    printf '%s\tskipped-already-done\n' "$name" >> "$MANIFEST"
done
for vcf in ${todo[@]+"${todo[@]}"}; do
    name="$(chunk_name "$vcf")"
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
    n_tsv=0
    hdr=""
    for tsv in "$OUTDIR"/*.acmg_plp.tsv; do
        [ "$tsv" = "$OUT_ALL" ] && continue
        [ -s "$tsv" ] || continue
        # Every chunk must agree on the column layout, or `tail -n +2` splices
        # rows under the wrong headers and the misalignment is invisible.
        this_hdr="$(head -1 "$tsv")"
        if [ "$first" = 1 ]; then
            hdr="$this_hdr"; cat "$tsv" >> "$OUT_ALL"; first=0
        else
            if [ "$this_hdr" != "$hdr" ]; then
                echo "[batch] ERROR: header mismatch in $tsv; refusing to concatenate" >&2
                echo "  expected: $hdr" >&2
                echo "  found:    $this_hdr" >&2
                exit 4
            fi
            tail -n +2 "$tsv" >> "$OUT_ALL"
        fi
        n_tsv=$((n_tsv+1))
    done
    # A concatenation over fewer chunks than were asked for is the dropped-chunk
    # incident again: the table is well-formed, the cohort is quietly smaller.
    if [ "$n_tsv" -ne "${#inputs[@]}" ]; then
        echo "[batch] ERROR: concatenated $n_tsv ACMG table(s) but $((${#inputs[@]})) input VCF(s) were given." >&2
        echo "        Missing chunks are listed as anything other than 'ok'/'skipped-already-done' in $MANIFEST." >&2
        echo "        Re-run the failed chunks first; $OUT_ALL would under-count carriers." >&2
        exit 4
    fi
    echo "[batch] concatenated ACMG table: $OUT_ALL ($n_tsv chunks)" >&2
fi

[ "$n_fail" -eq 0 ] || exit 1
