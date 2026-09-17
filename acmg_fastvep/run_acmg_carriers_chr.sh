#!/usr/bin/env bash
# Per-chromosome driver for build_acmg_carriers.py.
#
# Runs the retrieval-only MANE-select script over every *.fastvep.vcf.gz
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
USAGE
    exit 2
}

HERE="$(cd "$(dirname "$0")" && pwd)"
INDIR=""; OUTDIR=""; JOBS=4
while [ $# -gt 0 ]; do
    case "$1" in
        --in-dir)    INDIR="$2"; shift 2 ;;
        -o|--outdir) OUTDIR="$2"; shift 2 ;;
        -j|--jobs)   JOBS="$2"; shift 2 ;;
        -h|--help)   usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$INDIR" ] && [ -n "$OUTDIR" ] || usage
[ -d "$INDIR" ] || { echo "ERROR: --in-dir not a directory: $INDIR" >&2; exit 3; }

shopt -s nullglob
chunks=("$INDIR"/*.fastvep.vcf.gz)
shopt -u nullglob
[ ${#chunks[@]} -gt 0 ] || { echo "ERROR: no *.fastvep.vcf.gz in $INDIR" >&2; exit 3; }
echo "[acmg-carriers] ${#chunks[@]} chunk(s) in $INDIR" >&2

mkdir -p "$OUTDIR"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/parts"

cat > "$WORK/run_one.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
vcf="\$1"
[ -s "\$vcf" ] || { echo "ERROR: empty or missing chunk: \$vcf" >&2; exit 3; }
name="\$(basename "\$vcf")"; name="\${name%.fastvep.vcf.gz}"
python3 "$HERE/build_acmg_carriers.py" --vcf "\$vcf" \\
    --out-variants "$WORK/parts/\${name}.variants.tsv" \\
    --out-carriers "$WORK/parts/\${name}.carriers.tsv"
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
    exit 4
fi

echo "[acmg-carriers] done: $OUTDIR/acmg_variants.tsv ($(( $(wc -l < "$OUTDIR/acmg_variants.tsv") - 1 )) rows), " \
     "$OUTDIR/acmg_carriers.tsv ($(( $(wc -l < "$OUTDIR/acmg_carriers.tsv") - 1 )) rows)" >&2
