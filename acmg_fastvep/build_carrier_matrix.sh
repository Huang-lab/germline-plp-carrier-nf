#!/usr/bin/env bash
# Standalone carrier-matrix assembly (no VEP/fastVEP recompute).
#
# Combines, for ONE chromosome/run, the ClinVar + fastVEP-ACMG per-variant
# classifications with genotypes from the QC'd VCFs into a single per-person
# carrier matrix — reusing the pipeline's bin/build_carrier_matrix.py. This is
# the cheap way to fold standalone fastVEP ACMG results into the ClinVar carrier
# output without re-running the whole Nextflow pipeline.
#
# Inputs (per chromosome) — reuses files you already have:
#   --norm-dir     results-chrN/norm_qc            (*.norm.vcf.gz  → genotypes)
#   --clinvar-dir  results-chrN/variants           (*.clinvar_plp.tsv)
#   --acmg-dir     results-fastvep-chrN            (*.acmg_plp.tsv)
#
# Output: one carrier_matrix.tsv (long): chr,pos,ref,alt,gene,person_id,GT,
#   zygosity,is_clinvar_PLP,is_acmg_PLP,is_AM_PLP  — carriers of any P/LP variant
#   under either framework.
set -euo pipefail

usage() {
    cat >&2 <<'U'
Usage: build_carrier_matrix.sh --norm-dir <norm_qc/> -o <out.tsv> [options]

Required:
  --norm-dir <dir>     dir of QC'd *.norm.vcf.gz (genotype source)
  -o, --out <file>     output combined carrier_matrix.tsv
  (and at least one of:)
  --clinvar-dir <dir>  dir of *.clinvar_plp.tsv   (e.g. results-chrN/variants)
  --acmg-dir <dir>     dir of *.acmg_plp.tsv       (e.g. results-fastvep-chrN)

Options:
  --keep <file>        sample keep-list (one ID per line)
  --bcftools <bin>     bcftools binary (default: bcftools on PATH; else python fallback)
  --out-wide <file>    also write a wide pivot
  -j, --jobs <N>       parallel per-chunk GT extractions (default 4; set to the
                       number of cores the job reserves)
  --allow-partial      proceed even when a chunk has no classifier TSV (default:
                       refuse — a missing chunk silently under-counts carriers)
U
    exit 2
}

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
NORM=""; OUT=""; CVDIR=""; ACDIR=""; KEEP=""; BCF="bcftools"; WIDE=""; JOBS="${JOBS:-4}"
ALLOW_PARTIAL=0
while [ $# -gt 0 ]; do
    case "$1" in
        --norm-dir)     NORM="$2"; shift 2 ;;
        -o|--out)       OUT="$2"; shift 2 ;;
        --clinvar-dir)  CVDIR="$2"; shift 2 ;;
        --acmg-dir)     ACDIR="$2"; shift 2 ;;
        --keep)         KEEP="$2"; shift 2 ;;
        --bcftools)     BCF="$2"; shift 2 ;;
        --out-wide)     WIDE="$2"; shift 2 ;;
        -j|--jobs)      JOBS="$2"; shift 2 ;;
        --allow-partial) ALLOW_PARTIAL=1; shift ;;
        -h|--help)      usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$NORM" ] && [ -n "$OUT" ] || usage
[ -d "$NORM" ] || { echo "ERROR: --norm-dir not a directory: $NORM" >&2; exit 3; }
[ -n "$CVDIR$ACDIR" ] || { echo "ERROR: need --clinvar-dir and/or --acmg-dir" >&2; exit 3; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
mkdir -p "$(dirname "$OUT")"

# concat per-chunk TSVs into one, keeping a single header.
#
# Writes the number of chunks it consumed to $CONCAT_N so the caller can check
# it against the chunk count. Without that check this function is the standalone
# path's version of the batch3 dropped-chunk incident: the Nextflow pipeline
# joins classifier TSVs to their VCF on chunk_id, so a missing chunk cannot go
# unnoticed there, whereas here a glob over whatever happens to be on disk
# produces a well-formed table over a quietly smaller cohort.
CONCAT_N=0
concat() {   # dir  glob  out
    local dir="$1" g="$2" out="$3" first=1 hdr="" this_hdr=""
    CONCAT_N=0
    : > "$out"
    shopt -s nullglob
    for f in "$dir"/$g; do
        [ -s "$f" ] || continue
        this_hdr="$(head -1 "$f")"
        if [ "$first" = 1 ]; then
            hdr="$this_hdr"; cat "$f" >> "$out"; first=0
        else
            # `tail -n +2` assumes every chunk has the same columns in the same
            # order. If one does not, the rows land under the wrong headers and
            # nothing downstream can see it.
            if [ "$this_hdr" != "$hdr" ]; then
                echo "ERROR: header mismatch in $f" >&2
                echo "  expected: $hdr" >&2
                echo "  found:    $this_hdr" >&2
                return 1
            fi
            tail -n +2 "$f" >> "$out"
        fi
        CONCAT_N=$((CONCAT_N+1))
    done
    shopt -u nullglob
    [ -s "$out" ]
}

# Refuse to build a carrier matrix from fewer chunks than the cohort has, unless
# the caller says that is what they meant.
check_chunk_count() {   # label  n_found  n_expected  dir
    local label="$1" found="$2" expected="$3" dir="$4"
    [ "$found" -eq "$expected" ] && return 0
    if [ "$ALLOW_PARTIAL" = 1 ]; then
        echo "[carrier] WARNING: $label has $found chunk table(s) for $expected chunk VCF(s) in $dir." >&2
        echo "          --allow-partial given, continuing; carrier counts are for the chunks present only." >&2
        return 0
    fi
    echo "ERROR: $label has $found chunk table(s) but --norm-dir holds $expected chunk VCF(s)." >&2
    echo "       Missing: $(( expected - found )). A carrier matrix built now would be well-formed" >&2
    echo "       and would under-count carriers, with nothing downstream able to tell." >&2
    echo "       Finish or re-run the missing chunks in $dir, or pass --allow-partial if a" >&2
    echo "       subset is genuinely what you want." >&2
    exit 4
}

# Chunk inventory first: the number of QC'd VCFs is the cohort's chunk count,
# and every classifier directory is expected to hold one table per chunk.
shopt -s nullglob
vcfs=("$NORM"/*.norm.vcf.gz)
[ ${#vcfs[@]} -eq 0 ] && vcfs=("$NORM"/*.vcf.gz)
[ ${#vcfs[@]} -eq 0 ] && vcfs=("$NORM"/*.norm.vcf)
[ ${#vcfs[@]} -eq 0 ] && vcfs=("$NORM"/*.vcf)
shopt -u nullglob
[ ${#vcfs[@]} -gt 0 ] || { echo "ERROR: no VCFs in $NORM" >&2; exit 3; }
n_chunks=${#vcfs[@]}
echo "[carrier] $n_chunks chunk VCF(s) in $NORM" >&2

CV=""; AC=""
if [ -n "$CVDIR" ]; then
    concat "$CVDIR" "*.clinvar_plp.tsv" "$WORK/clinvar.tsv" && CV="$WORK/clinvar.tsv" || true
    check_chunk_count "--clinvar-dir" "$CONCAT_N" "$n_chunks" "$CVDIR"
fi
if [ -n "$ACDIR" ]; then
    concat "$ACDIR" "*.acmg_plp.tsv" "$WORK/acmg.tsv" && AC="$WORK/acmg.tsv" || true
    check_chunk_count "--acmg-dir" "$CONCAT_N" "$n_chunks" "$ACDIR"
fi
[ -n "$CV$AC" ] || { echo "ERROR: no classifier TSVs found under the given dirs" >&2; exit 3; }

# qualifying positions = union of is_clinvar_PLP==1 and is_acmg_PLP==1 (header-located columns)
: > "$WORK/qual.txt"
for spec in ${CV:+"$CV:is_clinvar_PLP"} ${AC:+"$AC:is_acmg_PLP"}; do
    f="${spec%%:*}"; flag="${spec##*:}"
    awk -F'\t' -v flag="$flag" '
        NR==1 { for (i=1;i<=NF;i++) c[$i]=i; next }
        (c[flag] && $(c[flag])==1) { print $(c["chr"])"\t"$(c["pos"]) }
    ' "$f" >> "$WORK/qual.txt"
done
sort -u "$WORK/qual.txt" > "$WORK/qual.pos.txt"
echo "[carrier] $(wc -l < "$WORK/qual.pos.txt") qualifying positions" >&2

# genotypes at qualifying positions from every norm VCF → gt.tsv (chr pos ref alt sample GT)
: > "$WORK/gt.tsv"

if [ -s "$WORK/qual.pos.txt" ]; then
    have_bcf=0; command -v "$BCF" >/dev/null 2>&1 && have_bcf=1
    if [ "$have_bcf" = 1 ]; then
        echo "[carrier] GT extraction via bcftools (C speed): ${#vcfs[@]} VCF(s)" >&2
    else
        echo "[carrier] WARNING: '$BCF' not found — using the pure-Python fallback," >&2
        echo "          which is 10-50x slower. Load bcftools (ml bcftools) for real data." >&2
    fi
    # Per-chunk extraction runs in PARALLEL (each chunk is an independent full
    # scan of its own VCF) — the dominant cost. Each writes its own part file;
    # order is irrelevant since build_carrier_matrix.py aggregates.
    mkdir -p "$WORK/parts"
    cat > "$WORK/extract_one.sh" <<EXTRACT
#!/usr/bin/env bash
set -euo pipefail
vcf="\$1"
part="$WORK/parts/\$(basename "\$vcf").gt"
if [ "$have_bcf" = 1 ]; then
    "$BCF" query -T "$WORK/qual.pos.txt" \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n' "\$vcf" \\
      | awk -F'\\t' 'BEGIN{OFS="\\t"}{for(i=5;i<=NF;i++){split(\$i,kv,"=");print \$1,\$2,\$3,\$4,kv[1],kv[2]}}' \\
      > "\$part"
else
    python3 "$REPO/bin/fixture_gt_extract.py" --vcf "\$vcf" --positions "$WORK/qual.pos.txt" --out "\$part"
fi
echo "[gt] done \$(basename "\$vcf")" >&2
EXTRACT
    chmod +x "$WORK/extract_one.sh"

    echo "[carrier] extracting ${#vcfs[@]} chunk(s) with ${JOBS} parallel worker(s)" >&2
    printf '%s\n' "${vcfs[@]}" | xargs -r -P "$JOBS" -I{} "$WORK/extract_one.sh" {}
    # One part file per chunk, or a chunk's genotypes are missing from the
    # matrix. `|| true` here previously turned that into a silently smaller
    # cohort; xargs reports a failed worker, but an absent part file after a
    # successful run would not have been noticed.
    shopt -s nullglob
    parts=("$WORK"/parts/*.gt)
    shopt -u nullglob
    if [ ${#parts[@]} -ne "$n_chunks" ]; then
        echo "ERROR: GT extraction produced ${#parts[@]} part file(s) for $n_chunks chunk(s)." >&2
        exit 4
    fi
    cat "${parts[@]}" >> "$WORK/gt.tsv"
fi

python3 "$REPO/bin/build_carrier_matrix.py" \
    --gt "$WORK/gt.tsv" \
    ${CV:+--clinvar "$CV"} ${AC:+--acmg "$AC"} \
    ${KEEP:+--keep "$KEEP"} \
    --out-long "$OUT" ${WIDE:+--out-wide "$WIDE"}

echo "[done] $OUT  ($(( $(wc -l < "$OUT") - 1 )) carrier rows)" >&2
