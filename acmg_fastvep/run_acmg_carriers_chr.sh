#!/usr/bin/env bash
# Per-chromosome driver for build_acmg_carriers.py / build_acmg_carriers_pick.py.
#
# Runs the chosen carriers script over every *.fastvep.vcf.gz
# chunk in one chromosome's directory IN PARALLEL (xargs -P — same pattern
# already used in run_fastvep_batch.sh's GT extraction and
# build_carrier_matrix.sh's per-chunk work), then concatenates the per-chunk
# outputs into the chromosome's two final tables. No fastVEP re-run: this
# only reads annotation fastVEP already computed.
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage: run_acmg_carriers_chr.sh --in-dir <results-fastvep-chrN_v6/> -o <outdir> [options]

Required:
  --in-dir <dir>    directory of *.fastvep.vcf.gz chunks for ONE chromosome
  -o, --outdir <dir> output directory (created if absent)

Options:
  -j, --jobs <N>    parallel chunk workers (default 4)
  --mode <fastvep_pick|mane|pick>
                    fastvep_pick (default; fastVEP v0.5.0 output): build_acmg_carriers_pick.py
                      --selection fastvep_pick, i.e. fastVEP's own PICK=1 block per
                      (allele, gene); fails on output without a PICK field.
                    mane: build_acmg_carriers.py, MANE Select only.
                    pick: DEPRECATED build_acmg_carriers_pick.py --selection python_pick
                      (Python port of the pick order), for comparing with v6 output.
  --pick-order <s>  (pick mode) passed through as --pick-order
  --keep-flanking-genes
                    (pick / fastvep_pick) keep upstream/downstream-only gene rows
  --allow-pick-anomalies
                    (fastvep_pick) exit 0 even if some (allele, gene) had no or several
                    PICK=1 blocks. Default: tables are written, then exit 5.
USAGE
    exit 2
}

HERE="$(cd "$(dirname "$0")" && pwd)"
ALLOW_ANOM=0; INDIR=""; OUTDIR=""; JOBS=4; MODE=fastvep_pick; PICK_ORDER=""; KEEP_FLANKING=0
while [ $# -gt 0 ]; do
    case "$1" in
        --in-dir)    INDIR="$2"; shift 2 ;;
        -o|--outdir) OUTDIR="$2"; shift 2 ;;
        -j|--jobs)   JOBS="$2"; shift 2 ;;
        --mode)      MODE="$2"; shift 2 ;;
        --pick-order) PICK_ORDER="$2"; shift 2 ;;
        --keep-flanking-genes) KEEP_FLANKING=1; shift ;;
        --allow-pick-anomalies) ALLOW_ANOM=1; shift ;;
        -h|--help)   usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$INDIR" ] && [ -n "$OUTDIR" ] || usage
[ -d "$INDIR" ] || { echo "ERROR: --in-dir not a directory: $INDIR" >&2; exit 3; }
# EXTRA: shell-quoted options pasted verbatim into the generated per-chunk
# script below (where \${name} is that chunk's name).
EXTRA=""
case "$MODE" in
    mane) SCRIPT="build_acmg_carriers.py"
          if [ -n "$PICK_ORDER" ] || [ "$KEEP_FLANKING" = 1 ]; then
              echo "ERROR: --pick-order/--keep-flanking-genes need --mode pick" >&2; exit 2
          fi ;;
    pick) SCRIPT="build_acmg_carriers_pick.py"
          echo "[acmg-carriers] WARN: --mode pick (python_pick) is DEPRECATED; fastVEP v0.5.0" \
               "output should use --mode fastvep_pick" >&2
          EXTRA='--selection python_pick --no-markers --stats-json "$PARTS/${name}.stats.json"'
          if [ -n "$PICK_ORDER" ]; then EXTRA="$EXTRA --pick-order $(printf '%q' "$PICK_ORDER")"; fi
          if [ "$KEEP_FLANKING" = 1 ]; then EXTRA="$EXTRA --keep-flanking-genes"; fi ;;
    fastvep_pick)
          SCRIPT="build_acmg_carriers_pick.py"
          if [ -n "$PICK_ORDER" ]; then
              echo "ERROR: --pick-order needs --mode pick; fastvep_pick uses fastVEP's own PICK flag" >&2; exit 2
          fi
          # Each chunk is let through; the chromosome total decides, after concat.
          EXTRA='--selection fastvep_pick --allow-pick-anomalies --no-markers --stats-json "$PARTS/${name}.stats.json"'
          if [ "$KEEP_FLANKING" = 1 ]; then EXTRA="$EXTRA --keep-flanking-genes"; fi ;;
    *) echo "ERROR: --mode must be fastvep_pick, mane or pick (got $MODE)" >&2; exit 2 ;;
esac
# Never write into the directory being read: the fastVEP outputs are inputs only.
# Checked before anything is created, so a bad --outdir leaves --in-dir untouched.
in_real="$(realpath -e "$INDIR")"
out_real="$(realpath -m "$OUTDIR")"
case "$out_real/" in
    "$in_real"/*) echo "ERROR: --outdir ($OUTDIR) is inside --in-dir ($INDIR); refusing" >&2; exit 3 ;;
esac
mkdir -p "$OUTDIR"

shopt -s nullglob
chunks=("$INDIR"/*.fastvep.vcf.gz)
shopt -u nullglob
[ ${#chunks[@]} -gt 0 ] || { echo "ERROR: no *.fastvep.vcf.gz in $INDIR" >&2; exit 3; }
echo "[acmg-carriers] mode=$MODE ${#chunks[@]} chunk(s) in $INDIR" >&2

WORK="$(mktemp -d)"
# Outcome marker in OUTDIR (carrier_markers.py): neither while running, FAILED
# (reason, counts, exit code, time) on any non-zero exit — the tables may be
# there and look complete — DONE on success. Readers refuse a FAILED dir.
python3 "$HERE/carrier_markers.py" clear "$OUTDIR"
FAIL_REASON="a step failed (see this job's stderr)"
on_exit() {
    local rc=$?
    rm -rf "$WORK"
    if [ "$rc" != 0 ]; then
        local counts=""
        [ -s "$OUTDIR/summary.txt" ] && counts="$(sed -n 's/^ *fastVEP PICK flag: //p' "$OUTDIR/summary.txt")"
        python3 "$HERE/carrier_markers.py" failed "$OUTDIR" "$rc" "$FAIL_REASON" "$counts" || true
        echo "[acmg-carriers] wrote $OUTDIR/FAILED (exit $rc: $FAIL_REASON)" >&2
    fi
}
trap on_exit EXIT
trap 'FAIL_REASON="killed by a signal"; exit 143' TERM INT
mkdir -p "$WORK/parts"

cat > "$WORK/run_one.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
vcf="\$1"
[ -s "\$vcf" ] || { echo "ERROR: empty or missing chunk: \$vcf" >&2; exit 3; }
name="\$(basename "\$vcf")"; name="\${name%.fastvep.vcf.gz}"
PARTS="$WORK/parts"
python3 "$HERE/$SCRIPT" --vcf "\$vcf" \\
    --out-variants "$WORK/parts/\${name}.variants.tsv" \\
    --out-carriers "$WORK/parts/\${name}.carriers.tsv" \\
    $EXTRA 2> "$WORK/parts/\${name}.stderr" \\
    || { cat "$WORK/parts/\${name}.stderr" >&2; exit 1; }
echo "[acmg-carriers] done \$(basename "\$vcf")" >&2
EOF
chmod +x "$WORK/run_one.sh"

printf '%s\n' "${chunks[@]}" | xargs -r -P "$JOBS" -I{} "$WORK/run_one.sh" {}

# Concatenate, refusing on header mismatch or a missing chunk output — the
# same safety pattern build_carrier_matrix.sh uses (the batch3 dropped-chunk
# incident: a well-formed table over a quietly smaller cohort is worse than
# an error).
concat() {   # glob  out
    local g="$1" out="$2" first=1 hdr="" this_hdr="" n=0
    : > "$out"
    shopt -s nullglob
    for f in "$WORK"/parts/$g; do
        [ -s "$f" ] || continue
        this_hdr="$(head -1 "$f")"
        if [ "$first" = 1 ]; then
            hdr="$this_hdr"; cat "$f" >> "$out"; first=0
        else
            if [ "$this_hdr" != "$hdr" ]; then
                echo "ERROR: header mismatch in $f" >&2
                echo "  expected: $hdr" >&2
                echo "  found:    $this_hdr" >&2
                return 1
            fi
            tail -n +2 "$f" >> "$out"
        fi
        n=$((n+1))
    done
    shopt -u nullglob
    echo "$n"
}

n_variants=$(concat "*.variants.tsv" "$OUTDIR/acmg_variants.tsv")
n_carriers=$(concat "*.carriers.tsv" "$OUTDIR/acmg_carriers.tsv")

if [ "$n_variants" -ne "${#chunks[@]}" ] || [ "$n_carriers" -ne "${#chunks[@]}" ]; then
    echo "ERROR: expected ${#chunks[@]} chunk output(s) but concatenated $n_variants variant" \
         "table(s) and $n_carriers carrier table(s). A chunk's output is missing — refusing" \
         "to publish a carrier table that would silently under-count the cohort." >&2
    FAIL_REASON="chunk outputs missing: $n_variants variant / $n_carriers carrier tables for ${#chunks[@]} chunks"
    exit 4
fi

if [ "$MODE" = pick ] || [ "$MODE" = fastvep_pick ]; then
    # One summary for the whole chromosome (sum of every chunk's counters),
    # plus each chunk's own stderr, kept next to the tables.
    agg_args=("$WORK"/parts/*.stats.json)
    [ "$ALLOW_ANOM" = 1 ] && agg_args+=(--allow-pick-anomalies)
    set +e
    python3 "$HERE/build_acmg_carriers_pick.py" --aggregate "${agg_args[@]}" > "$OUTDIR/summary.txt"
    agg_rc=$?
    set -e
    cat "$OUTDIR/summary.txt" >&2
    case "$agg_rc" in 0|3) ;; *) echo "ERROR: aggregating chunk stats failed ($agg_rc)" >&2
        FAIL_REASON="aggregating chunk stats failed ($agg_rc)"; exit 4 ;; esac
fi
cat "$WORK"/parts/*.stderr > "$OUTDIR/chunks.stderr.log"

if [ "${agg_rc:-0}" = 3 ]; then
    echo "ERROR: fastvep_pick found (allele, gene) rows without exactly one PICK=1 (see" \
         "summary.txt and chunks.stderr.log). Tables are written; failing so this is looked at;" \
         "--allow-pick-anomalies accepts it." >&2
    FAIL_REASON="PICK anomalies: (allele, gene) rows without exactly one PICK=1"
    exit 5
fi
python3 "$HERE/carrier_markers.py" "done" "$OUTDIR" \
    "mode=$MODE allow_pick_anomalies=$ALLOW_ANOM $(sed -n 's/^ *fastVEP PICK flag: //p' "$OUTDIR/summary.txt" 2>/dev/null || true)"
echo "[acmg-carriers] done: $OUTDIR/acmg_variants.tsv ($(( $(wc -l < "$OUTDIR/acmg_variants.tsv") - 1 )) rows), " \
     "$OUTDIR/acmg_carriers.tsv ($(( $(wc -l < "$OUTDIR/acmg_carriers.tsv") - 1 )) rows)" >&2
