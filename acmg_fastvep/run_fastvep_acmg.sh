#!/usr/bin/env bash
# Back-compat shim: ACMG-AMP classification via fastVEP.
# Equivalent to `run_fastvep.sh --acmg <args...>`. See run_fastvep.sh for the
# full option list. Kept so existing docs/commands that call this name still work.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
exec "$HERE/run_fastvep.sh" --acmg "$@"
