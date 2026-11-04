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
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: submit_acmg_carriers.sh --batch-root <dir> [options]

Required:
  --batch-root <dir>   dir containing results-fastvep-chr*_v6/ subdirs

Options:
  --only <chrspec>      only submit the chromosome dir(s) matching this glob
                        fragment (e.g. "chr21", "chr21_v6") — for a first test
                        run before submitting the full array
  --out-root <dir>      where results-carriers[-pick]-chrN_v6/ dirs, the array
                        manifest and logs/ are written (default: --batch-root
                        in mane mode; REQUIRED and must lie outside
                        --batch-root in pick mode, which never writes there)
  --mode <mane|pick>    mane (default): MANE Select only (build_acmg_carriers.py)
                        pick: per-gene pick fallback (build_acmg_carriers_pick.py)
  --pick-order <s>      pick mode: tier order (default in the script)
  --keep-flanking-genes pick mode: keep upstream/downstream-only gene rows
  -j, --jobs <N>        parallel chunk workers per chromosome task (default 4)
  -P, --project <alloc> LSF allocation (default: $MINERVA_ALLOCATION)
  -h, --help
USAGE
    exit 2
}

HERE="$(cd "$(dirname "$0")" && pwd)"
BATCH_ROOT=""; ONLY=""; OUT_ROOT=""; JOBS=4; MODE=mane; PICK_ORDER=""; KEEP_FLANKING=0
PROJECT="${MINERVA_ALLOCATION:-}"
while [ $# -gt 0 ]; do
    case "$1" in
        --batch-root) BATCH_ROOT="$2"; shift 2 ;;
        --only)       ONLY="$2"; shift 2 ;;
        --out-root)   OUT_ROOT="$2"; shift 2 ;;
        -j|--jobs)    JOBS="$2"; shift 2 ;;
        --mode)       MODE="$2"; shift 2 ;;
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
    pick) [ -n "$OUT_ROOT" ] || { echo "ERROR: --mode pick needs an explicit --out-root" >&2; exit 2; }
          b_real="$(realpath -e "$BATCH_ROOT")"; o_real="$(realpath -m "$OUT_ROOT")"
          case "$o_real/" in
              "$b_real"/*) echo "ERROR: --out-root ($OUT_ROOT) is inside --batch-root; the fastVEP" \
                                "outputs are read-only for this run — pick a separate dir" >&2; exit 3 ;;
          esac
          # Fail at submit time, not in 24 array tasks, on a bad order.
          python3 -c 'import sys; sys.path[:0]=[sys.argv[1], sys.argv[1]+"/.."]; import build_acmg_carriers_pick as P; P.parse_pick_order(sys.argv[2] or P.DEFAULT_PICK_ORDER)' \
              "$HERE" "$PICK_ORDER" || { echo "ERROR: bad --pick-order" >&2; exit 2; } ;;
    *) echo "ERROR: --mode must be mane or pick" >&2; exit 2 ;;
esac
mkdir -p "$OUT_ROOT"

# Discover chromosome directories. One glob, no per-chromosome loop of bsub
# calls — the manifest just records what the array will index into.
shopt -s nullglob
if [ -n "$ONLY" ]; then
    dirs=("$BATCH_ROOT"/results-fastvep-*"$ONLY"*/)
else
    dirs=("$BATCH_ROOT"/results-fastvep-chr*_v6/)
fi
shopt -u nullglob
[ ${#dirs[@]} -gt 0 ] || { echo "ERROR: no results-fastvep-chr*_v6/ dirs found under $BATCH_ROOT (--only=$ONLY)" >&2; exit 3; }

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
    bsub -P "$PROJECT" -J "acmg-carriers-$MODE[1-$N]" \
    -env "all, AC_MANIFEST=$MANIFEST, AC_OUT_ROOT=$OUT_ROOT, AC_JOBS=$JOBS, AC_REPO=$(cd "$HERE/.." && pwd)" \
    < "$HERE/run_acmg_carriers_task.lsf"

echo "[submit] submitted array acmg-carriers-$MODE[1-$N] for $N chromosome(s) -> $OUT_ROOT" >&2
echo "[submit] manifest kept at $MANIFEST (do not delete until the array finishes)." >&2
