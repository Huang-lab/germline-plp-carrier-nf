#!/usr/bin/env bash
# submit_fv050_pipeline.sh — submit the whole fastVEP v0.5.0 ACMG run in one go.
#
# For each chromosome it submits three LSF steps, chained so nothing waits on you:
#   1. annotate: one job per chunk  (run_fastvep_batch.sh --lsf; chunks already done
#                                    with the same provenance are skipped, not resubmitted)
#   2. verify:   -w "ended(<all chunk jobs of this chr>)"  re-runs the batch WITHOUT --lsf:
#                every finished chunk is skipped, a failed/killed one is redone here,
#                and the job exits non-zero if any chunk still fails
#   3. carriers: -w "done(verify)"  run_acmg_carriers_chr.sh --mode fastvep_pick
# A chromosome whose verify fails never builds carriers (its carrier job stays PEND;
# bkill it after fixing). Job names carry a run tag, so a re-submission never waits
# on an older run's jobs.
#
# Run from the login node, with the environment the jobs need already loaded:
#   module load python/3.12.5 htslib
#   export FV_ENV=$FV_HOME/env.sh
#   acmg_fastvep/submit_fv050_pipeline.sh --batch-dir $B --skip "13 17" --dry-run
#   acmg_fastvep/submit_fv050_pipeline.sh --batch-dir $B --skip "13 17" \
#       --fv-out-root /sc/arion/scratch/$USER/fv050
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: submit_fv050_pipeline.sh --batch-dir <B> [options]

  --batch-dir <B>       holds results-chr<N>/norm_qc/*.norm.vcf.gz (input)
  --chroms "<list>"     chromosomes (default: "1 2 ... 22 X Y"; missing ones are skipped)
  --skip "<list>"       chromosomes to leave out (e.g. "13 17", already done)
  --fv-out-root <dir>   where results-fastvep-chr<N>_fv050/ go (default: <B>).
                        Annotated VCFs are ~2.8 GB per chunk: use scratch if <B> is short.
  --carrier-root <dir>  where results-carriers-fastvep_pick-chr<N>_fv050/ go (default: <B>)
  --log-dir <dir>       LSF logs of the verify/carrier jobs (default: ./logs)
  --modules "<list>"    loaded inside verify/carrier jobs (default: "python/3.12.5 htslib")
  --run-tag <tag>       job-name tag (default: date +%m%d%H%M)
  --chunk-threads <N>   per chunk job (default 4)   --chunk-mem <MB> (default 16000)
  --chunk-walltime <hh:mm> (default 4:00)           --queue <q> (default premium)
  --dry-run             print every bsub instead of submitting

Environment: FV_ENV=<$FV_HOME/env.sh>, MINERVA_ALLOCATION.
USAGE
    exit 2
}

HERE="$(cd "$(dirname "$0")" && pwd)"
B=""; CHROMS="$(seq 1 22 | tr '\n' ' ')X Y"; SKIP=""; FVROOT=""; CARROOT=""
LOGDIR="$PWD/logs"; MODULES="python/3.12.5 htslib"; TAG="$(date +%m%d%H%M)"
CTHREADS=4; CMEM=16000; CWALL="4:00"; QUEUE="premium"; DRY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --batch-dir)      B="$2"; shift 2 ;;
        --chroms)         CHROMS="$2"; shift 2 ;;
        --skip)           SKIP="$2"; shift 2 ;;
        --fv-out-root)    FVROOT="$2"; shift 2 ;;
        --carrier-root)   CARROOT="$2"; shift 2 ;;
        --log-dir)        LOGDIR="$2"; shift 2 ;;
        --modules)        MODULES="$2"; shift 2 ;;
        --run-tag)        TAG="$2"; shift 2 ;;
        --chunk-threads)  CTHREADS="$2"; shift 2 ;;
        --chunk-mem)      CMEM="$2"; shift 2 ;;
        --chunk-walltime) CWALL="$2"; shift 2 ;;
        --queue)          QUEUE="$2"; shift 2 ;;
        --dry-run)        DRY=1; shift ;;
        -h|--help)        usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$B" ] || usage
[ -d "$B" ] || { echo "ERROR: --batch-dir not a directory: $B" >&2; exit 3; }
: "${FV_ENV:?export FV_ENV=\$FV_HOME/env.sh first}"
[ -s "$FV_ENV" ] || { echo "ERROR: FV_ENV not found: $FV_ENV" >&2; exit 3; }
# shellcheck disable=SC1090
. "$FV_ENV"
: "${MINERVA_ALLOCATION:?set MINERVA_ALLOCATION}"
FVROOT="${FVROOT:-$B}"; CARROOT="${CARROOT:-$B}"
mkdir -p "$LOGDIR"; LOGDIR="$(cd "$LOGDIR" && pwd)"

# The chunk jobs inherit this shell's environment (python3 parses the ACMG table in
# each one): refuse to submit 500 jobs that would all die on a broken python3.
python3 -c 'import sys; assert sys.version_info >= (3, 9)' 2>/dev/null || {
    echo "ERROR: python3 >= 3.9 does not run in this shell (module load python/3.12.5?)" >&2; exit 3; }
command -v tabix >/dev/null || { echo "ERROR: tabix not on PATH (module load htslib)" >&2; exit 3; }

if [ "$DRY" = 1 ]; then
    # Every bsub, including the ones run_fastvep_batch.sh makes, is printed instead.
    STUB="$(mktemp -d)"; trap 'rm -rf "$STUB"' EXIT
    printf '#!/usr/bin/env bash\nprintf "[dry-run] bsub"; printf " %%q" "$@"; echo\necho "Job <0> is submitted to queue <dry>."\n' > "$STUB/bsub"
    chmod +x "$STUB/bsub"; PATH="$STUB:$PATH"
fi

q() { printf '%q' "$1"; }
fv_args=(--fastvep "$FV_BIN" --gff3 "$FV_GFF3" --fasta "$FV_FASTA"
         --transcript-cache "$FV_CACHE" --hgvs --acmg --sa-dir "$FV_SA")
prelude="module load $MODULES && . $(q "$FV_ENV")"

# Jobs LSF knows about. A re-submission must not race an earlier run still in the
# queue (a duplicate chunk job hits the chunk lock, the new verify records it FAILED,
# and the new carrier job then pends forever), and a reused tag would make done(<verify>)
# also match the old run's verify job.
ACTIVE=""; KNOWN=""
if command -v bjobs >/dev/null 2>&1; then
    ACTIVE="$(bjobs -noheader -o job_name 2>/dev/null || true)"
    KNOWN="$(bjobs -a -noheader -o job_name 2>/dev/null || true)"
fi
if printf '%s\n' "$KNOWN" | grep -q "^fv050-$TAG-"; then
    echo "ERROR: LSF already has jobs named fv050-$TAG-*; pick another --run-tag" >&2; exit 3
fi

declare -a OK=() BAD=() MISSING=() BUSY=()
for c in $CHROMS; do
    case " $SKIP " in *" $c "*) echo "[pipeline] chr$c: skipped (--skip)" >&2; continue ;; esac
    IN="$B/results-chr$c/norm_qc"
    if ! ls "$IN"/*.norm.vcf.gz >/dev/null 2>&1; then
        echo "[pipeline] chr$c: no $IN/*.norm.vcf.gz, skipped" >&2; MISSING+=("$c"); continue
    fi
    if printf '%s\n' "$ACTIVE" | grep -Eq "^fv050-(.*-)?(acmg|verify|carriers)-chr${c}(-|\$)"; then
        echo "[pipeline] chr$c: jobs of an earlier run are still pending/running; not submitted" \
             "(bjobs | grep chr$c, or wait for them)" >&2
        BUSY+=("$c"); continue
    fi
    OUT="$FVROOT/results-fastvep-chr${c}_fv050"
    CAR="$CARROOT/results-carriers-fastvep_pick-chr${c}_fv050"
    P="fv050-$TAG-acmg-chr${c}-"
    VJ="fv050-$TAG-verify-chr$c"; CJ="fv050-$TAG-carriers-chr$c"
    slog="$LOGDIR/fv050-$TAG-submit-chr$c.log"

    # 1. chunk jobs
    if ! "$HERE/run_fastvep_batch.sh" --lsf -P "$MINERVA_ALLOCATION" --queue "$QUEUE" \
            --walltime "$CWALL" --mem "$CMEM" --threads "$CTHREADS" --job-prefix "$P" \
            --in-dir "$IN" -o "$OUT" "${fv_args[@]}" 2>&1 | tee "$slog" >&2; then
        echo "[pipeline] chr$c: chunk submission FAILED (see $slog); no verify/carrier job." \
             "Chunk jobs it did submit may still run: bkill -J '${P}*', or let them finish and" \
             "re-run this script for chr$c (finished chunks are skipped)." >&2
        BAD+=("$c"); continue
    fi
    n_sub="$(sed -n 's/.*submitted \([0-9][0-9]*\) LSF jobs.*/\1/p' "$slog" | tail -1)"
    n_sub="${n_sub:-0}"
    wait_chunks=(); [ "$n_sub" -gt 0 ] && wait_chunks=(-w "ended(${P}*)")

    # 2. verify (re-runs any chunk that did not finish; exits non-zero if one still fails)
    vcmd="$prelude && $(q "$HERE/run_fastvep_batch.sh") --threads 8 --in-dir $(q "$IN") -o $(q "$OUT")"
    for a in "${fv_args[@]}"; do vcmd+=" $(q "$a")"; done
    # rusage[mem] is per core on Minerva: 8 x 4000 = 32 GB, 4 x 4000 = 16 GB.
    if ! bsub -P "$MINERVA_ALLOCATION" -q "$QUEUE" -J "$VJ" -n 8 -R "rusage[mem=4000] span[hosts=1]" \
         -W 12:00 ${wait_chunks[@]+"${wait_chunks[@]}"} \
         -oo "$LOGDIR/$VJ.%J.out" -eo "$LOGDIR/$VJ.%J.err" "$vcmd"; then
        echo "[pipeline] chr$c: verify bsub FAILED; its $n_sub chunk job(s) run unverified." \
             "Re-run this script for chr$c once they end." >&2
        BAD+=("$c"); continue
    fi

    # 3. carriers, only after verify succeeded
    ccmd="$prelude && $(q "$HERE/run_acmg_carriers_chr.sh") --mode fastvep_pick -j 4"
    ccmd+=" --in-dir $(q "$OUT") -o $(q "$CAR")"
    if ! bsub -P "$MINERVA_ALLOCATION" -q "$QUEUE" -J "$CJ" -n 4 -R "rusage[mem=4000] span[hosts=1]" \
         -W 4:00 -w "done($VJ)" \
         -oo "$LOGDIR/$CJ.%J.out" -eo "$LOGDIR/$CJ.%J.err" "$ccmd"; then
        echo "[pipeline] chr$c: carriers bsub FAILED; after $VJ ends, run" \
             "run_acmg_carriers_chr.sh for chr$c by hand." >&2
        BAD+=("$c"); continue
    fi

    echo "[pipeline] chr$c: $n_sub chunk job(s) ${P}* -> $VJ -> $CJ" >&2
    OK+=("$c")
done

echo >&2
echo "[pipeline] run tag $TAG$([ "$DRY" = 1 ] && echo ' (DRY RUN: nothing submitted)')" >&2
echo "[pipeline] submitted: ${OK[*]:-none}" >&2
[ ${#MISSING[@]} -eq 0 ] || echo "[pipeline] no input:  ${MISSING[*]}" >&2
[ ${#BUSY[@]} -eq 0 ] || echo "[pipeline] busy (earlier jobs still queued), not submitted: ${BUSY[*]}" >&2
[ ${#BAD[@]} -eq 0 ] || { echo "[pipeline] FAILED to submit: ${BAD[*]}" >&2; exit 1; }
echo "[pipeline] watch: bjobs -J 'fv050-$TAG-*' | awk 'NR>1{print \$3}' | sort | uniq -c" >&2
