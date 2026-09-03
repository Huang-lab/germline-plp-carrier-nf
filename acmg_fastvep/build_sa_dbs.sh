#!/usr/bin/env bash
# Build fastVEP supplementary-annotation databases into $FV_ROOT/sa_db.
#
# Run on a COMPUTE node (offline — reads local source files in $FV_ROOT/sa_src).
# Whole-genome sources (clinvar, revel, gnomad_genes, clinvar_protein, omim,
# phylop, repeatmasker) are built once if missing. gnomAD and SpliceAI are
# per-chromosome / precomputed downloads; gnomAD auto-detects from the
# gnomad.exomes.*.sites.chrN.*.bgz files present in sa_src.
#
# v6 note: as of the v6 sa_db build (2026-08-19), this script only covered
# clinvar/revel/gnomad_genes/gnomad. The other 5 sources the ACMG classifier
# reads (spliceai, phylop, omim, clinvar_protein, repeatmasker) were built by
# hand directly against $FV_ROOT/sa_db and were not previously reflected here.
# This update adds them so the full 9-source stack is reproducible from this
# script instead of living only as untracked state on Minerva. See
# docs/ACMG_SETUP.md upstream (Huang-lab/fastVEP) for the full source list and
# what each one drives.
#
# omim: this repo builds real OMIM `genemap2.txt` (downloaded directly from
# omim.org, requires an OMIM license), not ClinGen Gene-Disease Validity. Both
# parse into the same .oga under the `omim` source key and are both supported
# by `fastvep sa-build --source omim`, but they are not equivalent: ClinGen GDV
# carries an explicit Refuted/Disputed tier that real OMIM's genemap2 does not,
# so a gene ClinGen has refuted may still pass the PVS1 gene-validity gate here
# where it would not against the upstream repo's own published ACMG benchmark
# (which defaults to ClinGen GDV). Worth knowing if a PVS1 call on a
# contested-validity gene ever looks surprising.
#
# Usage:
#   export FASTVEP=/.../fastvep/bin/fastvep
#   export FV_ROOT=/.../data/fastvep            # has sa_src/ ; sa_db/ is created
#   export CLINVAR=/.../clinvar_YYYYMMDD.vcf.gz  # optional (whole-genome clinvar source)
#   export CLINVAR_VARIANT_SUMMARY=/.../variant_summary.txt.gz  # optional (clinvar_protein);
#     defaults to the shared ClinVar refs copy at
#     /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/variants_PLP/data/refs/clinvar/variant_summary.txt.gz
#     if unset and nothing is found under sa_src/
#   export OMIM_GENEMAP2=/.../genemap2.txt       # optional (real OMIM; requires a license)
#   acmg_fastvep/build_sa_dbs.sh                 # build whole-genome + all gnomAD chrs found
#   acmg_fastvep/build_sa_dbs.sh 2 3 4           # build whole-genome (if missing) + gnomAD chr2,3,4
set -euo pipefail

FASTVEP="${FASTVEP:?set FASTVEP=<fastvep binary>}"
FV_ROOT="${FV_ROOT:?set FV_ROOT=<tool root with sa_src/>}"
SRC="$FV_ROOT/sa_src"
DB="$FV_ROOT/sa_db"
CLINVAR="${CLINVAR:-}"
CLINVAR_VARIANT_SUMMARY="${CLINVAR_VARIANT_SUMMARY:-}"
OMIM_GENEMAP2="${OMIM_GENEMAP2:-}"
mkdir -p "$DB"

command -v "$FASTVEP" >/dev/null 2>&1 || { echo "ERROR: fastvep not found: $FASTVEP" >&2; exit 3; }
[ -d "$SRC" ] || { echo "ERROR: sources dir missing: $SRC" >&2; exit 3; }

build_once() {   # source  outbase  inputfile
    local src="$1" out="$2" in="$3"
    if ls "$DB/$out".* >/dev/null 2>&1; then echo "[skip] $out already built"; return; fi
    [ -s "$in" ] || { echo "[skip] $out — no input file: $in" >&2; return; }
    echo "[build] $out  <-  $in"
    "$FASTVEP" sa-build --source "$src" -i "$in" -o "$DB/$out" --assembly GRCh38
}

# --- whole-genome sources (built once) ---
[ -n "$CLINVAR" ] && build_once clinvar clinvar "$CLINVAR"
build_once revel revel "$SRC/revel_with_transcript_ids"
GC="$(ls "$SRC"/gnomad*constraint*metrics*.tsv 2>/dev/null | head -1 || true)"
[ -n "$GC" ] && build_once gnomad_genes gnomad_genes "$GC"

# clinvar_protein: gene-level protein/splice index. Must be built from
# variant_summary.txt.gz, not the allele VCF — the VCF carries no HGVS c.
# token, so it yields few protein entries and no splice index at all.
# For the v6 build this lived outside sa_src, under the shared ClinVar refs
# dir rather than this tool's own source tree — default to that path, since
# nothing under sa_src/ will match the fallback glob.
VS="${CLINVAR_VARIANT_SUMMARY:-$(ls "$SRC"/variant_summary.txt.gz 2>/dev/null | head -1 || true)}"
[ -z "$VS" ] && [ -s "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/variants_PLP/data/refs/clinvar/variant_summary.txt.gz" ] && \
    VS="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/variants_PLP/data/refs/clinvar/variant_summary.txt.gz"
[ -n "$VS" ] && build_once clinvar_protein clinvar_protein "$VS"

# omim: real OMIM genemap2.txt (see the note above). Real OMIM requires a
# license; ClinGen Gene-Disease Validity is the free/default alternative
# upstream's own docs recommend if you don't have one.
OM="${OMIM_GENEMAP2:-$(ls "$SRC"/genemap2.txt 2>/dev/null | head -1 || true)}"
[ -n "$OM" ] && build_once omim omim "$OM"

# repeatmasker: interval-level (.osi), built from a BED. This repo's sa_src
# ships both a raw UCSC rmsk.txt.gz and repeatmasker_to_bed.py to convert it;
# run that conversion first if repeatmasker.bed isn't already present.
RM_BED="$(ls "$SRC"/repeatmasker.bed 2>/dev/null | head -1 || true)"
if [ -z "$RM_BED" ] && [ -s "$SRC/rmsk.txt.gz" ] && [ -x "$SRC/repeatmasker_to_bed.py" ]; then
    echo "[build] repeatmasker.bed  <-  $SRC/rmsk.txt.gz (via repeatmasker_to_bed.py)"
    python3 "$SRC/repeatmasker_to_bed.py" "$SRC/rmsk.txt.gz" > "$SRC/repeatmasker.bed"
    RM_BED="$SRC/repeatmasker.bed"
fi
[ -n "$RM_BED" ] && build_once custom_bed repeatmasker "$RM_BED"

# phylop: allele-level conservation. Upstream's recommended path distills this
# (and SpliceAI) directly out of gnomAD v4.1's own INFO columns via
# scripts/extract_gnomad_scores.py (Huang-lab/fastVEP), which needs no UCSC
# download. This repo's sa_src instead ships raw per-chromosome UCSC
# phyloP100way wigFix files (chrN.phyloP100way.wigFix.gz, plus a combined
# hg38.phyloP100way.wigFix.gz) — build from whichever fastVEP's phylop source
# parser accepts; adjust the input glob below if it expects the combined file
# rather than per-chromosome ones.
PHYLOP_SRC="$(ls "$SRC"/hg38.phyloP100way.wigFix.gz 2>/dev/null | head -1 || true)"
[ -n "$PHYLOP_SRC" ] && build_once phylop phylop "$PHYLOP_SRC"

# spliceai: allele-level, extremely dense (near every base of every gene body
# scores 3+ alternates). Built at fastVEP's default sa-build chunk width
# (--chunk-bits 20), a source this dense silently fails most lookups — see
# Huang-lab/fastVEP#101, fixed upstream but not yet in this repo's pinned
# fastvep binary as of the v6 build. Until that fix is deployed here, build at
# a narrower chunk width via the two-step sa-build(--format osa) + sa-convert
# path documented in fastVEP's docs/SUPPLEMENTARY_ANNOTATIONS.md, not a plain
# build_once call:
#   "$FASTVEP" sa-build --source spliceai -i "$SRC"/spliceai_scores.*.vcf.gz \
#       -o /tmp/spliceai_v1 --assembly GRCh38 --format osa
#   "$FASTVEP" sa-convert -i /tmp/spliceai_v1.osa -o "$DB/spliceai" --chunk-bits 16
SPLICEAI_SRC="$(ls "$SRC"/spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz 2>/dev/null | head -1 || true)"
if [ -n "$SPLICEAI_SRC" ] && ! ls "$DB"/spliceai.* >/dev/null 2>&1; then
    echo "[skip] spliceai — needs the sa-build(--format osa) + sa-convert(--chunk-bits) path above, not build_once; see comment"
fi

# --- gnomAD per-chromosome ---
chrs=("$@")
if [ ${#chrs[@]} -eq 0 ]; then
    shopt -s nullglob
    for f in "$SRC"/gnomad.exomes.*.sites.chr*.vcf.bgz; do
        c="$(basename "$f" | sed -E 's/.*\.sites\.chr([0-9XY]+)\.vcf\.bgz/\1/')"
        chrs+=("$c")
    done
    shopt -u nullglob
fi
for c in "${chrs[@]}"; do
    f="$(ls "$SRC"/gnomad.exomes.*.sites.chr${c}.vcf.bgz 2>/dev/null | head -1 || true)"
    [ -n "$f" ] || { echo "[skip] no gnomAD source for chr$c in $SRC" >&2; continue; }
    build_once gnomad "gnomad_chr${c}" "$f"
done

echo "[done] databases in $DB:"
ls -lh "$DB"
