#!/usr/bin/env bash
# Verify a fastVEP v0.5.0 ACMG sa_db before anything is annotated against it,
# and seal it with a manifest that later runs can check cheaply.
#
# fastVEP does not fail on an incomplete --sa-dir: it classifies without the
# missing evidence and says so once on stderr (docs/ACMG_SETUP.md). So "the
# files exist" is checked here as exactly as it can be, and then whether each
# source answers a query is checked by fastVEP's own scripts/check_acmg_stack.py.
#
# Full mode (default; run as a job, check_acmg_stack.py annotates with the GFF3,
# ~16 GB):
#   1. the file set is EXACTLY the expected 100 databases, each non-empty with a
#      .done marker from build_sa_dbs.sh; any other .osa/.osa2/.osa.idx/.oga/.osi
#      in sa_db is a FAIL (fastVEP loads every one it finds: a stray
#      gnomad_genomes_chr1.osa2 would silently replace exome frequencies);
#   2. every gnomad_chr* was built from gnomAD v4.1 EXOMES (never joint/genomes);
#   3. python3 $FV_SRC/scripts/check_acmg_stack.py --gff3 --fasta --sa-dir
#      --fastvep $FV_BIN --keep-json: each of the nine sources answers >= 1 probe;
#   4. omim.oga (ClinGen by default) knows a few established genes: by name in
#      the ClinGen TSV it was built from, and by the probe JSON's `omim` key;
#   5. writes sa_db/SA_DB_SHA256.txt (sha256sum format) + SA_DB_MANIFEST.tsv,
#      only after 1-4 pass.
# --quick (what run_fastvep*.sh call before a batch): step 1, plus each file's
#   size equal to the manifest's and none modified after the manifest. No
#   rehash, no annotation. Prints `sa_db_sha256 <hex of SA_DB_SHA256.txt>`.
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: check_sa_db.sh [--quick] [--sa-dir <dir>] [--skip-stack]

Environment (from `source $FV_HOME/env.sh`, or FV_ENV=<path to env.sh>):
  FV_SA (sa_db), and for the full check FV_BIN FV_GFF3 FV_FASTA FV_SRC FV_HOME

Options:
  --quick         file set + manifest check only (no hashing, no fastVEP run)
  --sa-dir <dir>  check this directory instead of $FV_SA
  --skip-stack    full mode without check_acmg_stack.py (e.g. no compute node yet);
                  the manifest is then NOT written, since the stack is unproven
  -h, --help
USAGE
    exit 2
}

if [ -n "${FV_ENV:-}" ]; then
    # shellcheck disable=SC1090
    . "$FV_ENV"
fi
QUICK=0; SKIP_STACK=0; DB="${FV_SA:-}"
while [ $# -gt 0 ]; do
    case "$1" in
        --quick)      QUICK=1; shift ;;
        --sa-dir)     DB="$2"; shift 2 ;;
        --skip-stack) SKIP_STACK=1; shift ;;
        -h|--help)    usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$DB" ] || { echo "ERROR: set FV_SA or pass --sa-dir" >&2; exit 2; }
[ -d "$DB" ] || { echo "ERROR: sa_db not a directory: $DB" >&2; exit 3; }

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}
size_of() { wc -c < "$1" | tr -d ' '; }

CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
EXPECTED="clinvar.osa2 clinvar_protein.oga gnomad_genes.oga omim.oga repeatmasker.osi"
for c in $CHROMS; do
    EXPECTED="$EXPECTED gnomad_chr${c}.osa2 spliceai_chr${c}.osa2 phylop_chr${c}.osa2"
    [ "$c" = Y ] || EXPECTED="$EXPECTED revel_chr${c}.osa2"   # REVEL has no chrY
done

FAIL=0
fail() { echo "FAIL: $*" >&2; FAIL=1; }

# --- 1. exact file set ---
n_exp=0
for f in $EXPECTED; do
    n_exp=$((n_exp+1))
    if [ ! -s "$DB/$f" ]; then fail "missing or empty: $f"
    elif [ ! -e "$DB/$f.done" ]; then fail "no .done marker (build did not finish?): $f"
    fi
done
shopt -s nullglob
for p in "$DB"/*.osa "$DB"/*.osa2 "$DB"/*.osa.idx "$DB"/*.oga "$DB"/*.osi; do
    b="$(basename "$p")"
    case " $EXPECTED " in *" $b "*) ;; *) fail "unexpected database file in sa_db (fastVEP would load it): $b" ;; esac
done
for d in "$DB"/*.done; do
    [ -e "${d%.done}" ] || fail "orphan marker without its database: $(basename "$d")"
done
shopt -u nullglob

# --- 2. gnomAD exomes only ---
for c in $CHROMS; do
    d="$DB/gnomad_chr${c}.osa2.done"
    [ -e "$d" ] || continue
    grep -q '^input .*gnomad\.exomes\.v4\.1\.sites\.chr' "$d" || \
        fail "gnomad_chr${c}.osa2 was not built from gnomAD v4.1 exomes ($(sed -n 's/^input  *//p' "$d"))"
done

[ "$FAIL" = 0 ] || { echo "check_sa_db: file-set check FAILED for $DB ($n_exp expected)" >&2; exit 1; }
echo "[check] file set OK: exactly $n_exp databases, each with .done" >&2

SHA="$DB/SA_DB_SHA256.txt"; MAN="$DB/SA_DB_MANIFEST.tsv"

if [ "$QUICK" = 1 ]; then
    [ -s "$SHA" ] && [ -s "$MAN" ] || { echo "FAIL: no $SHA / $MAN — run check_sa_db.sh (full) once" >&2; exit 1; }
    while IFS=$'\t' read -r f size _sha _rest; do
        [ "$f" = file ] && continue
        [ -s "$DB/$f" ] || { fail "manifest lists $f but it is gone"; continue; }
        [ "$(size_of "$DB/$f")" = "$size" ] || fail "$f size changed since the manifest"
        [ "$DB/$f" -nt "$SHA" ] && fail "$f modified after the manifest was written"
    done < "$MAN"
    n_man="$(($(wc -l < "$MAN") - 1))"
    [ "$n_man" = "$n_exp" ] || fail "manifest lists $n_man files, expected $n_exp"
    [ "$FAIL" = 0 ] || { echo "check_sa_db --quick: FAILED; re-run the full check_sa_db.sh" >&2; exit 1; }
    echo "sa_db_sha256 $(sha256_of "$SHA")"
    exit 0
fi

# --- 3. does each source answer? (fastVEP's own probe) ---
: "${FV_BIN:?set FV_BIN}"; : "${FV_SRC:?set FV_SRC}"; : "${FV_GFF3:?set FV_GFF3}"; : "${FV_FASTA:?set FV_FASTA}"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
if [ "$SKIP_STACK" = 1 ]; then
    echo "[check] --skip-stack: check_acmg_stack.py NOT run; manifest will not be written" >&2
else
    echo "[check] check_acmg_stack.py against $DB (annotates 8 ClinVar probes) ..." >&2
    if ! python3 "$FV_SRC/scripts/check_acmg_stack.py" --gff3 "$FV_GFF3" --fasta "$FV_FASTA" \
            --sa-dir "$DB" --fastvep "$FV_BIN" --keep-json "$WORK/probe.json" >&2; then
        fail "check_acmg_stack.py: a source in the stack answered nothing (see above)"
    fi
fi

# --- 4. omim.oga content ---
KNOWN_GENES="BRCA1 BRCA2 TP53 PTEN MYH7"
omim_in="$(sed -n 's/^input  *//p' "$DB/omim.oga.done")"
case "$omim_in" in
    *clingen_gdv.tsv)
        if [ -s "$omim_in" ]; then
            for g in $KNOWN_GENES; do
                awk -F'\t' -v g="$g" '$9 == g { found = 1 } END { exit !found }' "$omim_in" || \
                    fail "ClinGen TSV $omim_in (omim.oga's input) has no curation for $g"
            done
        else
            fail "omim.oga's input $omim_in is gone; cannot check its genes"
        fi ;;
    *) echo "[check] omim.oga built from $omim_in (not ClinGen); gene-name check skipped" >&2 ;;
esac
if [ -s "$WORK/probe.json" ]; then
    python3 - "$WORK/probe.json" <<'PY' || fail "omim key missing on probe genes BRCA1/PTEN/MYH7 (see above)"
import json, sys
variants = json.load(open(sys.argv[1]))
want = {"BRCA1_missense", "BRCA1_nonsense", "PTEN_missense", "MYH7_missense"}
missing = []
for v in variants:
    if v.get("id") in want and not any("omim" in ga for ga in (v.get("genes") or {}).values()):
        missing.append(v.get("id"))
seen = {v.get("id") for v in variants} & want
if missing or not seen:
    print(f"omim not attached for probes: {missing or sorted(want)}", file=sys.stderr)
    sys.exit(1)
print(f"[check] omim attached on probes: {sorted(seen)}", file=sys.stderr)
PY
fi

[ "$FAIL" = 0 ] || { echo "check_sa_db: FAILED for $DB — manifest not written" >&2; exit 1; }
if [ "$SKIP_STACK" = 1 ]; then
    echo "check_sa_db: file set + gene checks OK; stack NOT probed (--skip-stack); no manifest written" >&2
    exit 0
fi

# --- 5. seal ---
echo "[check] hashing $n_exp databases ..." >&2
: > "$WORK/sha"; printf 'file\tsize\tsha256\tbuilt_utc\tinput\tinput_sha256\n' > "$WORK/man"
# shellcheck disable=SC2086  # EXPECTED is a word list by design
for f in $(printf '%s\n' $EXPECTED | LC_ALL=C sort); do
    h="$(sha256_of "$DB/$f")"
    printf '%s  %s\n' "$h" "$f" >> "$WORK/sha"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$f" "$(size_of "$DB/$f")" "$h" \
        "$(sed -n 's/^built_utc  *//p' "$DB/$f.done")" "$(sed -n 's/^input  *//p' "$DB/$f.done")" \
        "$(sed -n 's/^input_sha256  *//p' "$DB/$f.done")" >> "$WORK/man"
done
mv -f "$WORK/man" "$MAN.partial"; mv -f "$WORK/sha" "$SHA.partial"
mv -f "$MAN.partial" "$MAN"; mv -f "$SHA.partial" "$SHA"
echo "[check] wrote $SHA and $MAN" >&2
echo "sa_db_sha256 $(sha256_of "$SHA")"
