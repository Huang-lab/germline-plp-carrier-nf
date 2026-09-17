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
  --out-root <dir>      where results-carriers-chrN_v6/ dirs are written
                        (default: same as --batch-root)
  -j, --jobs <N>        parallel chunk workers per chromosome task (default 4)
  -P, --project <alloc> LSF allocation (default: $MINERVA_ALLOCATION)
  -h, --help
USAGE
    exit 2
}

HERE="$(cd "$(dirname "$0")" && pwd)"
BATCH_ROOT=""; ONLY=""; OUT_ROOT=""; JOBS=4
PROJECT="${MINERVA_ALLOCATION:-}"
while [ $# -gt 0 ]; do
    case "$1" in
        --batch-root) BATCH_ROOT="$2"; shift 2 ;;
        --only)       ONLY="$2"; shift 2 ;;
        --out-root)   OUT_ROOT="$2"; shift 2 ;;
        -j|--jobs)    JOBS="$2"; shift 2 ;;
        -P|--project) PROJECT="$2"; shift 2 ;;
        -h|--help)    usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$BATCH_ROOT" ] || usage
[ -d "$BATCH_ROOT" ] || { echo "ERROR: --batch-root not a directory: $BATCH_ROOT" >&2; exit 3; }
[ -n "$PROJECT" ] || { echo "ERROR: need -P/--project or \$MINERVA_ALLOCATION" >&2; exit 3; }
OUT_ROOT="${OUT_ROOT:-$BATCH_ROOT}"

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

MANIFEST="$BATCH_ROOT/.acmg_carriers_manifest.$$.txt"
: > "$MANIFEST"
for d in "${dirs[@]}"; do
    printf '%s\n' "${d%/}" >> "$MANIFEST"
done
N=${#dirs[@]}
echo "[submit] $N chromosome dir(s) -> $MANIFEST" >&2
cat "$MANIFEST" >&2

mkdir -p "$BATCH_ROOT/logs"
AC_MANIFEST="$MANIFEST" AC_OUT_ROOT="$OUT_ROOT" AC_JOBS="$JOBS" AC_REPO="$(cd "$HERE/.." && pwd)" \
    bsub -P "$PROJECT" -J "acmg-carriers[1-$N]" \
    -env "all, AC_MANIFEST=$MANIFEST, AC_OUT_ROOT=$OUT_ROOT, AC_JOBS=$JOBS, AC_REPO=$(cd "$HERE/.." && pwd)" \
    < "$HERE/run_acmg_carriers_task.lsf"

echo "[submit] submitted array acmg-carriers[1-$N] for $N chromosome(s)." >&2
echo "[submit] manifest kept at $MANIFEST (do not delete until the array finishes)." >&2
