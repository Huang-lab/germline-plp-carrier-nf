#!/usr/bin/env bash
# Set up ANNOVAR + InterVar on Minerva.
#
# ANNOVAR is registration-gated. You must supply the tarball path yourself;
# this script never fetches it. InterVar is freely available (MIT).
#
# The ClinVar humandb version installed here MUST match the ClinVar VCF release
# fetched by setup/fetch_references.sh (pass the same $CLINVAR_RELEASE).

set -euo pipefail

: "${ANNOVAR_TARBALL:?Set ANNOVAR_TARBALL to the path of your registration-gated ANNOVAR tarball}"
: "${RESOURCES_DIR:?Set RESOURCES_DIR (same as for fetch_references.sh)}"
: "${CLINVAR_RELEASE:?Set CLINVAR_RELEASE — must match the ClinVar VCF release used by fetch_references.sh}"

if [ ! -s "$ANNOVAR_TARBALL" ]; then
    echo "ERROR: ANNOVAR tarball not found at $ANNOVAR_TARBALL" >&2
    echo "  Obtain it from https://annovar.openbioinformatics.org/ (registration required)." >&2
    exit 2
fi

ANNOVAR_DIR="$RESOURCES_DIR/annovar"
HUMANDB="$RESOURCES_DIR/humandb"
INTERVAR_DIR="$RESOURCES_DIR/InterVar"

mkdir -p "$RESOURCES_DIR" "$HUMANDB"

if [ ! -x "$ANNOVAR_DIR/annotate_variation.pl" ]; then
    tar -xzf "$ANNOVAR_TARBALL" -C "$RESOURCES_DIR"
fi

# Download ANNOVAR humandb DBs InterVar needs.
# NOTE: exact DB names on Minerva evolve; verify the current names with
#   `$ANNOVAR_DIR/annotate_variation.pl -downdb -webfrom annovar -h`
# and update below if newer releases have appeared.
DBS=( refGene "clinvar_${CLINVAR_RELEASE}" gnomad41_exome dbnsfp47a avsnp150 )
for db in "${DBS[@]}"; do
    if [ ! -e "$HUMANDB/hg38_${db}.txt" ] && [ ! -e "$HUMANDB/hg38_${db}.txt.gz" ]; then
        echo "downdb: $db" >&2
        "$ANNOVAR_DIR/annotate_variation.pl" -buildver hg38 -downdb -webfrom annovar "$db" "$HUMANDB"
    fi
done

# Clone InterVar (freely available).
if [ ! -x "$INTERVAR_DIR/Intervar.py" ]; then
    git clone --depth=1 https://github.com/WGLab/InterVar.git "$INTERVAR_DIR"
fi

# Write config.ini pointing InterVar at the ANNOVAR install + our humandb + intervardb.
cat > "$INTERVAR_DIR/config.ini" <<CFG
[InterVar]
buildver = hg38
inputfile_type = AVinput
current_version = 2.2.2
convert2annovar = $ANNOVAR_DIR/convert2annovar.pl
table_annovar = $ANNOVAR_DIR/table_annovar.pl
annotate_variation = $ANNOVAR_DIR/annotate_variation.pl
database_locat = $HUMANDB
database_intervar = $INTERVAR_DIR/intervardb
disease_variant_flag = clinvar_${CLINVAR_RELEASE}
CFG

echo "annovar_intervar=setup-complete clinvar=${CLINVAR_RELEASE}" >> "$RESOURCES_DIR/versions.txt"
echo "Done. Point params/<cohort>.yaml at:" >&2
echo "  annovar_dir      = $ANNOVAR_DIR" >&2
echo "  annovar_humandb  = $HUMANDB" >&2
echo "  intervar_dir     = $INTERVAR_DIR" >&2
echo "  intervar_config  = $INTERVAR_DIR/config.ini" >&2
echo "REMINDER: supply the AlphaMissense gene-specific calibration table (Chen/Pejaver 2026) yourself" >&2
echo "  and point params.am_calibration_tsv at it." >&2
