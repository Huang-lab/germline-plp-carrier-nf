#!/usr/bin/env bash
# Add gnomAD SpliceAI indel scores to the SNV-only SpliceAI database.
#
# The v6 spliceai.osa2 was built from Illumina's
# spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz alone, so indels
# get no SpliceAI score. Illumina's indel file is behind a BaseSpace account.
# gnomAD carries `spliceai_ds_max` on its own records, indels included
# (chr22 of v4.0 exomes: 299,149 indel scores of 3,527,809), so this takes the
# indels from there and merges them into the SNV file before one build.
#
# What the gnomAD half carries (see Huang-lab/fastVEP scripts/extract_gnomad_scores.py):
# the max delta score in all four DS fields and 0 in the DP position fields;
# fastVEP's ACMG criteria read only the max. Only indels gnomAD has seen get a
# score. SNVs keep Illumina's full record.
#
# Modes (everything under $WORK; the live $FV_ROOT/sa_db is never written):
#   extract <chr>  tabix-index the gnomAD chr if needed, extract SpliceAI, keep
#                  indels, drop the `chr` prefix to match the Illumina file
#   merge          merge-sort SNVs + all indels into one VCF (order checked),
#                  then sa-build(--format osa) + sa-convert(--chunk-bits 16) as
#                  in build_sa_dbs.sh, into $OUT/spliceai.*
#   check          on chr22 SNVs present in both sources, compare gnomAD's
#                  spliceai_ds_max with Illumina's max(DS): run before swapping
#   submit         LSF array job for the 24 extracts, merge after it ends
#   submit-merge   resubmit only the merge (e.g. after rerunning one extract)
#
# Usage (Minerva login node, with htslib and python3 loaded):
#   export MINERVA_ALLOCATION=acc_<project>
#   export FV_ROOT=/sc/arion/projects/.../data/fastvep   # has bin/ sa_src/ tools/
#   acmg_fastvep/build_spliceai_indels.sh submit
#   acmg_fastvep/build_spliceai_indels.sh check          # after extract 22 is done
# Then compare $OUT/spliceai.* against $FV_ROOT/sa_db/spliceai.* and swap by hand.
# `check` reads up to chr22 of the unindexed 27 GB SNV file unless it has a .tbi:
# run it as a job (bsub ... bash <this script> check). If sa-build/sa-convert
# fails, clear $OUT before rerunning `merge`.
#
# Optional overrides: FASTVEP, EXTRACT (extract_gnomad_scores.py), PYTHON (an
# interpreter that starts on compute nodes, e.g. /usr/bin/python3), GNOMAD_FMT
# (printf format of the per-chr gnomAD file name), WORK, OUT, MERGE_MEM (MB per
# core for the merge job; size it from the original spliceai sa-build's log).
set -euo pipefail
export LC_ALL=C

FV_ROOT="${FV_ROOT:?set FV_ROOT}"
FASTVEP="${FASTVEP:-$FV_ROOT/bin/fastvep}"
EXTRACT="${EXTRACT:-$FV_ROOT/tools/extract_gnomad_scores.py}"
PYTHON="${PYTHON:-python3}"
SRC="$FV_ROOT/sa_src"
SNV="$SRC/spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz"
GNOMAD_FMT="${GNOMAD_FMT:-gnomad.exomes.v4.0.sites.chr%s.vcf.bgz}"
WORK="${WORK:-$FV_ROOT/spliceai_indel_work}"
OUT="${OUT:-$FV_ROOT/sa_db_spliceai_indel}"
MERGE_MEM="${MERGE_MEM:-8000}"
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

# Skip headers and rank a chromosome the way the Illumina file orders them: 1..22, X, Y.
RANK_AWK='BEGIN{FS=OFS="\t"} /^#/{next} {c=$1; r=(c=="X")?23:(c=="Y")?24:(c ~ /^[0-9]+$/)?c+0:99; print r, $0}'

die() { echo "$*" >&2; exit 1; }

extract() {
    local c="$1" vcf
    vcf="$SRC/$(printf "$GNOMAD_FMT" "$c")"
    [ -f "$vcf" ] || die "missing $vcf"
    # A job killed mid-index can leave a partial .tbi; only trust one that lists the contig.
    if ! grep -qx "chr$c" <(tabix -l "$vcf" 2>/dev/null); then
        rm -f "$vcf.tbi" "$vcf.csi"
        tabix -p vcf "$vcf"
    fi
    mkdir -p "$WORK/chr$c"
    "$PYTHON" "$EXTRACT" --chrom "$c" --url "$vcf" --out-dir "$WORK/chr$c"
    awk 'BEGIN{FS=OFS="\t"} !/^#/ && length($4)!=length($5) {sub(/^chr/,"",$1); print}' \
        "$WORK/chr$c/spliceai_chr$c.vcf" > "$WORK/indels_chr$c.tmp"
    mv "$WORK/indels_chr$c.tmp" "$WORK/indels_chr$c.tsv"
    # chr22's full extract is kept for `check`
    if [ "$c" = 22 ]; then rm -f "$WORK/chr$c/phylop_chr$c.tsv"; else rm -rf "$WORK/chr$c"; fi
    echo "chr$c: $(wc -l < "$WORK/indels_chr$c.tsv") indels"
}

merge() {
    local c files=() combined="$WORK/spliceai_snv_plus_gnomad_indel.vcf.gz"
    local snv_fifo="$WORK/snv.fifo" indel_fifo="$WORK/indel.fifo" snv_pid indel_pid
    local live out_real
    live="$(CDPATH= cd "$FV_ROOT/sa_db" 2>/dev/null && pwd -P || true)"
    out_real="$(CDPATH= cd "$OUT" 2>/dev/null && pwd -P || true)"
    [ -z "$live" ] || [ "$live" != "$out_real" ] || die "OUT is the live sa_db: $OUT"
    ! ls "$OUT"/spliceai.* >/dev/null 2>&1 || die "$OUT already has spliceai.*; move it away first"
    for c in "${CHROMS[@]}"; do
        [ -f "$WORK/indels_chr$c.tsv" ] || die "missing indels_chr$c.tsv; rerun: $0 extract $c"
        files+=("$WORK/indels_chr$c.tsv")
    done
    [ -f "$SNV" ] || die "missing $SNV"
    command -v bgzip >/dev/null || die "need bgzip on PATH (module load htslib)"

    # Paths are expanded now: the trap fires after merge() returns and its locals are gone.
    # It also kills the producers: if the pipeline fails before sort opens a FIFO,
    # its writer would otherwise block in open() forever.
    trap "kill \$(jobs -p) 2>/dev/null || true; rm -f $(printf '%q ' "$combined.tmp" "$WORK/merge_counts.tmp" "$snv_fifo" "$indel_fifo")" EXIT
    rm -f "$snv_fifo" "$indel_fifo"; mkfifo "$snv_fifo" "$indel_fifo"
    # Producers run as background jobs so their exit status can be checked: a
    # truncated read inside <(...) would just end the input early, unnoticed.
    ( set -o pipefail; zcat "$SNV" | awk "$RANK_AWK" ) > "$snv_fifo" & snv_pid=$!
    ( set -o pipefail; cat "${files[@]}" | awk "$RANK_AWK" ) > "$indel_fifo" & indel_pid=$!

    echo "[merge] SNVs + gnomAD indels -> $combined"
    {
        { zcat "$SNV" || true; } | awk '/^##/{print; next} /^#CHROM/{
            print "##SpliceAI_indels=gnomAD v4.0 exomes spliceai_ds_max (extract_gnomad_scores.py): DS fields carry the max, DP fields are 0"
            print; exit}'
        sort -m -t $'\t' -k1,1n -k3,3n "$snv_fifo" "$indel_fifo" \
          | awk -v counts="$WORK/merge_counts.tmp" 'BEGIN{FS=OFS="\t"}
              # input is rank \t record; check order, count, strip the rank
              { if ($1 < pr || ($1 == pr && $3+0 < pp)) { print "out of order at " $2 ":" $3 > "/dev/stderr"; bad=1; exit 1 }
                pr=$1; pp=$3+0; n++; if (length($5)!=length($6)) ni++
                sub(/^[^\t]*\t/, ""); print }
              END { if (!bad) print "records\t" n "\nindels\t" ni > counts }'
    } | bgzip -@ "${LSB_DJOB_NUMPROC:-4}" > "$combined.tmp"
    wait "$snv_pid" || die "reading $SNV failed; merged file discarded"
    wait "$indel_pid" || die "reading the indel files failed; merged file discarded"
    mv "$combined.tmp" "$combined"
    mv "$WORK/merge_counts.tmp" "$WORK/merge_counts.txt"
    cat "$WORK/merge_counts.txt"

    echo "[build] spliceai (two-step, --chunk-bits 16) -> $OUT"
    mkdir -p "$OUT"
    "$FASTVEP" sa-build --source spliceai -i "$combined" -o "$WORK/spliceai_v1" --assembly GRCh38 --format osa
    "$FASTVEP" sa-convert -i "$WORK/spliceai_v1.osa" -o "$OUT/spliceai" --chunk-bits 16
    rm -f "$WORK"/spliceai_v1.osa*
    ( cd "$OUT" && sha256sum spliceai.* | tee "$OUT/SPLICEAI_SHA256.txt" )
    echo "[done] run '$0 check', then compare with $FV_ROOT/sa_db/spliceai.* before swapping"
}

# Are gnomAD's scores on Illumina's scale? Compare max(DS) on chr22 SNVs found in both.
check() {
    local g="$WORK/chr22/spliceai_chr22.vcf"
    [ -f "$g" ] || die "missing $g; run: $0 extract 22"
    # Illumina chr22 records: straight from the index if there is one, else scan up to chr22.
    if [ -f "$SNV.tbi" ] || [ -f "$SNV.csi" ]; then tabix "$SNV" 22
    else { zcat "$SNV" || true; } | awk 'BEGIN{FS="\t"} $1=="22"{print; seen=1; next} seen{exit}'
    fi | awk 'BEGIN{FS="\t"; OFS="\t"}
          NR==FNR { if (length($4)==length($5)) { split($8, f, "|"); g[$2":"$4":"$5]=f[3]+0 }; next }
          # Illumina can list one variant once per overlapping gene: keep its max over all of them.
          { k=$2":"$4":"$5; if (!(k in g)) next
            split($8, f, "|"); for (i=3;i<=6;i++) if (!(k in m) || f[i]+0 > m[k]) m[k]=f[i]+0 }
          END { for (k in m) { d=g[k]-m[k]; if (d<0) d=-d; n++; s+=d; if (d>0.015) far++
                               if (n<=5) print "  " k, "gnomad=" g[k], "illumina=" m[k] }
                if (!n) { print "no shared chr22 SNVs found"; exit 1 }
                printf "shared chr22 SNVs: %d  mean |diff|: %.4f  |diff|>0.015: %d (%.2f%%)\n", n, s/n, far+0, 100*far/n }' \
          <(grep -v '^#' "$g") -
}

submit() {
    local alloc="${MINERVA_ALLOCATION:?set MINERVA_ALLOCATION}" self jid
    self="$(CDPATH= cd "$(dirname "$0")" && pwd)/$(basename "$0")"
    command -v tabix >/dev/null && command -v bgzip >/dev/null && command -v sha256sum >/dev/null \
        || die "need tabix, bgzip, sha256sum on PATH (module load htslib)"
    # A module/conda python can start here and still fail on a compute node (libpython not found).
    "$PYTHON" -c 'import sys; assert sys.version_info >= (3, 6)' || die "$PYTHON does not run; set PYTHON=/usr/bin/python3"
    [ -f "$EXTRACT" ] || die "missing $EXTRACT"
    ! ls "$WORK"/indels_chr*.tsv >/dev/null 2>&1 || die "$WORK has indel files from an earlier run; move them away or use submit-merge"
    mkdir -p "$WORK/logs"
    jid=$(bsub -P "$alloc" -q premium -J "spliceai-indel[1-${#CHROMS[@]}]" -n 1 \
         -R "rusage[mem=8000]" -W 12:00 \
         -oo "$WORK/logs/extract.%I.out" -eo "$WORK/logs/extract.%I.err" \
         bash "$self" extract-index | sed -n 's/^Job <\([0-9]*\)>.*/\1/p')
    [ -n "$jid" ] || die "extract array submission failed"
    echo "extract array: job $jid"
    submit_merge "$jid"
}

# ended(), not done(): merge itself refuses to run if any chromosome's indel file is missing.
submit_merge() {
    local alloc="${MINERVA_ALLOCATION:?set MINERVA_ALLOCATION}" self dep=()
    self="$(CDPATH= cd "$(dirname "$0")" && pwd)/$(basename "$0")"
    [ -n "${1:-}" ] && dep=(-w "ended($1)")
    mkdir -p "$WORK/logs"
    bsub -P "$alloc" -q premium -J spliceai-indel-merge ${dep[@]+"${dep[@]}"} \
         -n 4 -R "rusage[mem=$MERGE_MEM] span[hosts=1]" -W 24:00 \
         -oo "$WORK/logs/merge.out" -eo "$WORK/logs/merge.err" \
         bash "$self" merge
}

case "${1:-}" in
    extract)       extract "${2:?chromosome, e.g. 22 or X}" ;;
    extract-index) extract "${CHROMS[$((LSB_JOBINDEX - 1))]}" ;;
    merge)         merge ;;
    check)         check ;;
    submit)        submit ;;
    submit-merge)  submit_merge ;;
    *) sed -n '2,36p' "$0"; exit 2 ;;
esac
