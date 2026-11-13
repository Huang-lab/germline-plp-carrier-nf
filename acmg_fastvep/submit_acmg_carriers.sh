#!/usr/bin/env bash
# Submit ONE LSF array job covering every chromosome under a batch root, e.g.
#   /sc/arion/projects/rg_huangk06/variants_PLP_MSM/batch_002/
# which holds results-fastvep-chr1_v6/ ... results-fastvep-chrY_v6/.
#
# This issues exactly one `bsub -J "...[1-N]"` call — LSF itself fans the N
# chromosomes out to N array tasks. There is no shell loop over chromosomes
# submitting one bsub each.
#
# Usage:
#   export MINERVA_ALLOCATION=acc_<project>
#   acmg_fastvep/submit_acmg_carriers.sh --batch-root <dir> [--only chr21]
#   acmg_fastvep/submit_acmg_carriers.sh --batch-root <dir> --mode pick \
#       --out-root <new dir, e.g. results-acmg-pick/> [--only chr21]
#   acmg_fastvep/submit_acmg_carriers.sh --batch-root <dir> --mode fastvep_pick \
#       --dir-suffix _fv050 --out-root <new dir> [--only chr21]   # fastVEP v0.5.0 output
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: submit_acmg_carriers.sh --batch-root <dir> [options]

Required:
  --batch-root <dir>   dir containing results-fastvep-chr*<suffix>/ subdirs

Options:
  --only <chrspec>      only submit this chromosome. An exact name ("chr1" or
                        "chr1_v6" -> results-fastvep-chr1_v6/) wins; otherwise
                        it is a glob fragment, as before. So "chr1" no longer
                        also matches chr10..chr19.
  --exclude <chrspec>   leave this chromosome out (exact name, as above);
                        repeatable. E.g. run the array --exclude chr1 and chr1
                        on its own with a longer --wall-time.
  --out-root <dir>      where results-carriers[-pick]-chrN_v6/ dirs, the array
                        manifest and logs/ are written (default: --batch-root
                        in mane mode; REQUIRED in pick mode: a folder that is
                        not --batch-root itself nor inside any of its
                        results-fastvep-*/results-carriers-* dirs, e.g.
                        <batch-root>/results-acmg-pick)
  --dir-suffix <s>      suffix of the results-fastvep-chrN<s>/ dirs (default _v6;
                        _fv050 for the fastVEP v0.5.0 runs)
  --mode <mane|pick|fastvep_pick>
                        mane (default): MANE Select only (build_acmg_carriers.py)
                        pick: DEPRECATED Python pick port (build_acmg_carriers_pick.py)
                        fastvep_pick: fastVEP's own PICK=1 block (v0.5.0 output);
                          like pick, needs --out-root
  --pick-order <s>      pick mode: tier order (default in the script)
  --keep-flanking-genes pick mode: keep upstream/downstream-only gene rows
  -n, --cores <N>       cores per chromosome task (default 4, from the task file)
  -j, --jobs <N>        parallel chunk workers per task (default: --cores)
  -W, --wall-time <H:MM> LSF run limit per task (default 4:00, from the task
                        file). chr1 (48 chunks) needs far more than chr21 (6).
  -P, --project <alloc> LSF allocation (default: $MINERVA_ALLOCATION)
  -h, --help
USAGE
    exit 2
}

HERE="$(cd "$(dirname "$0")" && pwd)"
SUFFIX="_v6"; BATCH_ROOT=""; ONLY=""; OUT_ROOT=""; JOBS=""; CORES=""; WALL=""; EXCLUDE=(); MODE=mane; PICK_ORDER=""; KEEP_FLANKING=0
PROJECT="${MINERVA_ALLOCATION:-}"
while [ $# -gt 0 ]; do
    case "$1" in
        --batch-root) BATCH_ROOT="$2"; shift 2 ;;
        --only)       ONLY="$2"; shift 2 ;;
        --out-root)   OUT_ROOT="$2"; shift 2 ;;
        -j|--jobs)    JOBS="$2"; shift 2 ;;
        -n|--cores)   CORES="$2"; shift 2 ;;
        -W|--wall-time) WALL="$2"; shift 2 ;;
        --exclude)    EXCLUDE+=("$2"); shift 2 ;;
        --mode)       MODE="$2"; shift 2 ;;
        --dir-suffix) SUFFIX="$2"; shift 2 ;;
        --pick-order) PICK_ORDER="$2"; shift 2 ;;
        --keep-flanking-genes) KEEP_FLANKING=1; shift ;;
        -P|--project) PROJECT="$2"; shift 2 ;;
        -h|--help)    usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$BATCH_ROOT" ] || usage
[ -d "$BATCH_ROOT" ] || { echo "ERROR: --batch-root not a directory: $BATCH_ROOT" >&2; exit 3; }
[ -n "$PROJECT" ] || { echo "ERROR: need -P/--project or \$MINERVA_ALLOCATION" >&2; exit 3; }
case "$MODE" in
    mane) if [ -n "$PICK_ORDER" ] || [ "$KEEP_FLANKING" = 1 ]; then
              echo "ERROR: --pick-order/--keep-flanking-genes need --mode pick" >&2; exit 2
          fi
          OUT_ROOT="${OUT_ROOT:-$BATCH_ROOT}" ;;
    pick|fastvep_pick)
          [ -n "$OUT_ROOT" ] || { echo "ERROR: --mode $MODE needs an explicit --out-root" >&2; exit 2; }
          if [ "$MODE" = fastvep_pick ] && [ -n "$PICK_ORDER" ]; then
              echo "ERROR: --pick-order needs --mode pick" >&2; exit 2
          fi
          # The fastVEP chunk dirs and the MANE-only result dirs are read-only
          # for a pick run. A separate subfolder of the batch root (e.g.
          # <batch>/results-acmg-pick) is fine; the batch root itself, or
          # anything inside those dirs, is not.
          b_real="$(realpath -e "$BATCH_ROOT")"; o_real="$(realpath -m "$OUT_ROOT")"
          if [ "$o_real" = "$b_real" ]; then
              echo "ERROR: --out-root is --batch-root itself; use a subfolder such as" \
                   "$BATCH_ROOT/results-acmg-pick" >&2; exit 3
          fi
          for prot in "$b_real"/results-fastvep-* "$b_real"/results-carriers-*; do
              [ -d "$prot" ] || continue
              case "$o_real/" in
                  "$prot"/*) echo "ERROR: --out-root ($OUT_ROOT) is inside $prot, which this run" \
                                  "must not write to" >&2; exit 3 ;;
              esac
          done
          # Fail at submit time, not in 24 array tasks, on a bad order.
          python3 -c 'import sys; sys.path[:0]=[sys.argv[1], sys.argv[1]+"/.."]; import build_acmg_carriers_pick as P; P.parse_pick_order(sys.argv[2] or P.DEFAULT_PICK_ORDER)' \
              "$HERE" "$PICK_ORDER" || { echo "ERROR: bad --pick-order" >&2; exit 2; } ;;
    *) echo "ERROR: --mode must be mane, pick or fastvep_pick" >&2; exit 2 ;;
esac
mkdir -p "$OUT_ROOT"

# Discover chromosome directories. One glob, no per-chromosome loop of bsub
# calls — the manifest just records what the array will index into.
# exact_dir <spec>: results-fastvep-<spec>/ or results-fastvep-<spec><suffix>/, if it exists.
exact_dir() {
    local c
    for c in "$BATCH_ROOT/results-fastvep-$1" "$BATCH_ROOT/results-fastvep-${1}${SUFFIX}"; do
        if [ -d "$c" ]; then printf '%s/\n' "$c"; return 0; fi
    done
    return 1
}
shopt -s nullglob
if [ -n "$ONLY" ]; then
    if d="$(exact_dir "$ONLY")"; then
        dirs=("$d")
    else
        dirs=("$BATCH_ROOT"/results-fastvep-*"$ONLY"*/)
    fi
else
    dirs=("$BATCH_ROOT"/results-fastvep-chr*"$SUFFIX"/)
fi
shopt -u nullglob
for x in ${EXCLUDE[@]+"${EXCLUDE[@]}"}; do
    xd="$(exact_dir "$x")" || { echo "ERROR: --exclude $x: no results-fastvep-${x}[${SUFFIX}]/ under $BATCH_ROOT" >&2; exit 3; }
    kept=()
    for d in "${dirs[@]}"; do [ "${d%/}" = "${xd%/}" ] || kept+=("$d"); done
    dirs=(${kept[@]+"${kept[@]}"})
done
[ ${#dirs[@]} -gt 0 ] || { echo "ERROR: no results-fastvep-chr*$SUFFIX/ dirs left under $BATCH_ROOT (--only=$ONLY --exclude=${EXCLUDE[*]-})" >&2; exit 3; }
JOBS="${JOBS:-${CORES:-4}}"
LIMITS=()
if [ -n "$CORES" ]; then LIMITS+=(-n "$CORES"); fi
if [ -n "$WALL" ]; then LIMITS+=(-W "$WALL"); fi

MANIFEST="$OUT_ROOT/.acmg_carriers_manifest.$MODE.$$.txt"
: > "$MANIFEST"
for d in "${dirs[@]}"; do
    printf '%s\n' "${d%/}" >> "$MANIFEST"
done
N=${#dirs[@]}
echo "[submit] $N chromosome dir(s) -> $MANIFEST" >&2
cat "$MANIFEST" >&2

mkdir -p logs   # the task's #BSUB -oo/-eo paths are relative to this (submit) dir
# AC_PICK_ORDER holds commas, which `bsub -env "..."` would split on, so the
# mode options travel only via the environment that "-env all" propagates.
AC_MANIFEST="$MANIFEST" AC_OUT_ROOT="$OUT_ROOT" AC_JOBS="$JOBS" AC_REPO="$(cd "$HERE/.." && pwd)" \
AC_MODE="$MODE" AC_PICK_ORDER="$PICK_ORDER" AC_KEEP_FLANKING="$KEEP_FLANKING" \
    bsub -P "$PROJECT" -J "acmg-carriers-$MODE[1-$N]" ${LIMITS[@]+"${LIMITS[@]}"} \
    -env "all, AC_MANIFEST=$MANIFEST, AC_OUT_ROOT=$OUT_ROOT, AC_JOBS=$JOBS, AC_REPO=$(cd "$HERE/.." && pwd)" \
    < "$HERE/run_acmg_carriers_task.lsf"

echo "[submit] submitted array acmg-carriers-$MODE[1-$N] for $N chromosome(s) -> $OUT_ROOT" \
     "(cores=${CORES:-4 (task default)} jobs=$JOBS wall=${WALL:-4:00 (task default)})" >&2
echo "[submit] manifest kept at $MANIFEST (do not delete until the array finishes)." >&2
