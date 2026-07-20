#!/usr/bin/env bash
# Fetch freely-downloadable references for germline-plp-carrier-nf.
#
# INTENDED HOST: Minerva (behind the site HTTPS proxy).
# NOT INTENDED to be run in the authoring environment.
#
# Reuses existing lab paths when present; only downloads what's missing.
# Writes all outputs under $RESOURCES_DIR (a gitignored directory).
# Records versions to $RESOURCES_DIR/versions.txt for manifest.json provenance.

set -euo pipefail

: "${RESOURCES_DIR:?Set RESOURCES_DIR to a writable dir on Minerva scratch/work}"
: "${CLINVAR_RELEASE:?Set CLINVAR_RELEASE to a dated release, e.g. 20260401 — must match ANNOVAR humandb clinvar_<date>}"
: "${VEP_CACHE_VERSION:=113}"
: "${GNOMAD_VERSION:=v4.1}"

# Optional: reuse an existing lab-wide reference tree. Anything already at
# these paths is used as-is (no re-download).
: "${SHARED_REFS:=/sc/arion/projects/YOUR_PROJECT/refs}"

# Per-classifier gating — set to 1 to skip the corresponding resource block.
# For a ClinVar-only run, set:  SKIP_AM=1 SKIP_LOFTEE=1 SKIP_GNOMAD=1
: "${SKIP_AM:=0}"
: "${SKIP_LOFTEE:=0}"
: "${SKIP_GNOMAD:=0}"

mkdir -p "$RESOURCES_DIR"/{clinvar,vep_cache,vep_fasta,alphamissense,loftee,gnomad,plugins}
VERSIONS="$RESOURCES_DIR/versions.txt"
: >"$VERSIONS"

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*" >&2; }
reuse_or_download() {
    # $1 destination path, $2 shared candidate path, $3 URL
    local dest="$1" shared="$2" url="$3"
    if [ -s "$dest" ]; then log "keep existing $dest"; return; fi
    if [ -n "$shared" ] && [ -s "$shared" ]; then
        log "symlink shared → $dest"
        ln -sfn "$shared" "$dest"
        return
    fi
    log "download $url → $dest"
    curl -fsSL --retry 3 --retry-delay 5 -o "$dest.part" "$url"
    mv "$dest.part" "$dest"
}

# --- ClinVar (pinned date; MUST match ANNOVAR humandb clinvar_<date>) ---
CV_BASE="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/${CLINVAR_RELEASE:0:4}"
reuse_or_download "$RESOURCES_DIR/clinvar/clinvar_${CLINVAR_RELEASE}.vcf.gz" \
    "$SHARED_REFS/clinvar/clinvar_${CLINVAR_RELEASE}.vcf.gz" \
    "${CV_BASE}/clinvar_${CLINVAR_RELEASE}.vcf.gz"
reuse_or_download "$RESOURCES_DIR/clinvar/clinvar_${CLINVAR_RELEASE}.vcf.gz.tbi" \
    "$SHARED_REFS/clinvar/clinvar_${CLINVAR_RELEASE}.vcf.gz.tbi" \
    "${CV_BASE}/clinvar_${CLINVAR_RELEASE}.vcf.gz.tbi"
echo "clinvar=${CLINVAR_RELEASE}" >> "$VERSIONS"

# --- VEP cache + FASTA (huge; strongly prefer shared) ---
VEP_CACHE_TARBALL="homo_sapiens_vep_${VEP_CACHE_VERSION}_GRCh38.tar.gz"
reuse_or_download "$RESOURCES_DIR/vep_cache/${VEP_CACHE_TARBALL}" \
    "$SHARED_REFS/vep_cache/${VEP_CACHE_TARBALL}" \
    "https://ftp.ensembl.org/pub/release-${VEP_CACHE_VERSION}/variation/indexed_vep_cache/${VEP_CACHE_TARBALL}"
if [ ! -d "$RESOURCES_DIR/vep_cache/homo_sapiens" ]; then
    tar -xzf "$RESOURCES_DIR/vep_cache/${VEP_CACHE_TARBALL}" -C "$RESOURCES_DIR/vep_cache/"
fi
FASTA_NAME="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
reuse_or_download "$RESOURCES_DIR/vep_fasta/${FASTA_NAME}" \
    "$SHARED_REFS/vep_fasta/${FASTA_NAME}" \
    "https://ftp.ensembl.org/pub/release-${VEP_CACHE_VERSION}/fasta/homo_sapiens/dna/${FASTA_NAME}"
echo "vep_cache=${VEP_CACHE_VERSION}" >> "$VERSIONS"

# --- AlphaMissense data ---
if [ "$SKIP_AM" = "1" ]; then
    log "SKIP_AM=1 — skipping AlphaMissense downloads"
else
    AM_URL="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
    reuse_or_download "$RESOURCES_DIR/alphamissense/AlphaMissense_hg38.tsv.gz" \
        "$SHARED_REFS/alphamissense/AlphaMissense_hg38.tsv.gz" \
        "$AM_URL"
    # AlphaMissense does NOT publish a .tbi — must be generated with tabix.
    # We defer that to the annotate.sif container:
    #   singularity exec annotate.sif tabix -s 1 -b 2 -e 2 -S 1 \
    #       $RESOURCES_DIR/alphamissense/AlphaMissense_hg38.tsv.gz
    if [ -s "$SHARED_REFS/alphamissense/AlphaMissense_hg38.tsv.gz.tbi" ]; then
        ln -sfn "$SHARED_REFS/alphamissense/AlphaMissense_hg38.tsv.gz.tbi" \
                "$RESOURCES_DIR/alphamissense/AlphaMissense_hg38.tsv.gz.tbi"
    else
        log "note: AlphaMissense_hg38.tsv.gz.tbi not present; generate it later with:"
        log "  singularity exec annotate.sif tabix -s 1 -b 2 -e 2 -S 1 \\"
        log "      $RESOURCES_DIR/alphamissense/AlphaMissense_hg38.tsv.gz"
    fi
    reuse_or_download "$RESOURCES_DIR/plugins/AlphaMissense.pm" "" \
        "https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/${VEP_CACHE_VERSION}/AlphaMissense.pm"
    echo "alphamissense=hg38-2023" >> "$VERSIONS"
fi

# --- LOFTEE data (GRCh38) ---
if [ "$SKIP_LOFTEE" = "1" ]; then
    log "SKIP_LOFTEE=1 — skipping LOFTEE downloads"
else
    for f in human_ancestor.fa.gz human_ancestor.fa.gz.fai human_ancestor.fa.gz.gzi loftee.sql; do
        reuse_or_download "$RESOURCES_DIR/loftee/$f" \
            "$SHARED_REFS/loftee/$f" \
            "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/$f"
    done
    reuse_or_download "$RESOURCES_DIR/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw" \
        "$SHARED_REFS/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw" \
        "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
    echo "loftee=GRCh38" >> "$VERSIONS"
fi

# --- gnomAD v4 exome sites (~200GB; NEVER download if lab copy exists) ---
if [ "$SKIP_GNOMAD" = "1" ]; then
    log "SKIP_GNOMAD=1 — skipping gnomAD"
else
    GNOMAD_NAME="gnomad.exomes.${GNOMAD_VERSION}.sites.vcf.bgz"
    if [ -s "$SHARED_REFS/gnomad/${GNOMAD_NAME}" ]; then
        ln -sfn "$SHARED_REFS/gnomad/${GNOMAD_NAME}"     "$RESOURCES_DIR/gnomad/${GNOMAD_NAME}"
        ln -sfn "$SHARED_REFS/gnomad/${GNOMAD_NAME}.tbi" "$RESOURCES_DIR/gnomad/${GNOMAD_NAME}.tbi"
        log "reusing shared gnomAD at $SHARED_REFS/gnomad/${GNOMAD_NAME}"
    elif [ ! -s "$RESOURCES_DIR/gnomad/${GNOMAD_NAME}" ]; then
        log "WARNING: gnomAD not found locally or in \$SHARED_REFS."
        log "Set GNOMAD_URL_PREFIX to a mirror URL and re-run this script to fetch."
        if [ -n "${GNOMAD_URL_PREFIX:-}" ]; then
            curl -fsSL "$GNOMAD_URL_PREFIX/${GNOMAD_NAME}"     -o "$RESOURCES_DIR/gnomad/${GNOMAD_NAME}"
            curl -fsSL "$GNOMAD_URL_PREFIX/${GNOMAD_NAME}.tbi" -o "$RESOURCES_DIR/gnomad/${GNOMAD_NAME}.tbi"
        fi
    fi
    echo "gnomad=${GNOMAD_VERSION}" >> "$VERSIONS"
fi

log "versions:"
cat "$VERSIONS" >&2
log "done; write params/<cohort>.yaml paths to match $RESOURCES_DIR/*."
