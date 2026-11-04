#!/usr/bin/env bash
# One-time setup for the standalone fastVEP tool (Huang-lab/fastVEP).
#
# Builds/installs the `fastvep` binary and fetches the Ensembl GRCh38 GFF3 gene
# model needed for annotation. Reuses the pipeline's existing GRCh38 FASTA — it
# does NOT re-download the FASTA. Idempotent: steps whose outputs already exist
# are skipped. Nothing here submits jobs or writes into the repo.
#
# ACMG note: plain annotation needs only the GFF3 produced here. ACMG (--acmg)
# additionally needs supplementary databases built with `fastvep sa-build`; that
# heavy step is NOT automated here — see fastVEP docs/ACMG_SETUP.md.
#
# Version pinning (see Huang-lab/germline-plp-carrier-nf#3): a build with no
# --fastvep-ref used to clone and build whatever commit happened to be at the
# tip of Huang-lab/fastVEP's default branch that day, so two runs months apart
# could silently classify the same variant differently. There are no tagged
# releases upstream yet, so the pin is a commit SHA rather than a version tag.
# DEFAULT_FASTVEP_REF below is the fastVEP commit this script pins to unless
# overridden; bump it deliberately, not by accident, and note why in the commit
# message that changes it (e.g. "picks up #101, the SpliceAI dense-chunk fix").
DEFAULT_FASTVEP_REF="8c6b179f196a6dbfd60eb3773c1d7d69950db355"
set -euo pipefail

usage() {
    cat >&2 <<USAGE
Usage: setup_fastvep.sh --refdir <dir> [options]

Required:
  --refdir <dir>          where to put the Ensembl GFF3 (and its fastvep cache)

Options:
  --fasta <fa>            existing GRCh38 FASTA to reuse (from the pipeline).
                          If given and its .fai is missing, samtools faidx is run.
  --ensembl-release <N>   Ensembl release for the GFF3 (default: 115)
  --fastvep-src <dir>     path to a fastVEP checkout to build from.
                          If absent, the repo is git-cloned into <refdir>/fastVEP.
  --fastvep-ref <ref>     git ref (commit SHA, branch, or tag) of Huang-lab/fastVEP
                          to build. Default: $DEFAULT_FASTVEP_REF
                          Ignored with --skip-build, since no build happens.
  --skip-build            do not build fastvep (assume it is already on PATH)
  --build-cache           pre-build the transcript cache (else auto-built on 1st run)
  -h, --help
USAGE
    exit 2
}

REFDIR=""; FASTA=""; REL="115"; SRC=""; SKIP_BUILD=0; BUILD_CACHE=0
FASTVEP_REF="$DEFAULT_FASTVEP_REF"
while [ $# -gt 0 ]; do
    case "$1" in
        --refdir)          REFDIR="$2"; shift 2 ;;
        --fasta)           FASTA="$2"; shift 2 ;;
        --ensembl-release) REL="$2"; shift 2 ;;
        --fastvep-src)     SRC="$2"; shift 2 ;;
        --fastvep-ref)     FASTVEP_REF="$2"; shift 2 ;;
        --skip-build)      SKIP_BUILD=1; shift ;;
        --build-cache)     BUILD_CACHE=1; shift ;;
        -h|--help)         usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$REFDIR" ] || usage
mkdir -p "$REFDIR"

# --- 1. fastvep binary ---
BUILT_FRESH=0
if [ "$SKIP_BUILD" = 1 ]; then
    command -v fastvep >/dev/null 2>&1 || { echo "ERROR: --skip-build but fastvep not on PATH" >&2; exit 3; }
    echo "[setup] using existing fastvep: $(command -v fastvep)" >&2
    echo "[setup] WARN: --skip-build means this run trusts whatever fastvep is already" >&2
    echo "         on PATH — its provenance (git ref, build date) is unknown to this" >&2
    echo "         script and will not be recorded in $REFDIR/fastvep.provenance.txt" >&2
elif command -v fastvep >/dev/null 2>&1; then
    echo "[setup] fastvep already on PATH: $(command -v fastvep) ($(fastvep --version 2>/dev/null || echo '?'))" >&2
    echo "[setup] WARN: not rebuilt against --fastvep-ref=$FASTVEP_REF — an existing" >&2
    echo "         fastvep on PATH is used as-is. Remove it or pass --skip-build" >&2
    echo "         explicitly if that ambiguity matters for this run." >&2
else
    command -v cargo >/dev/null 2>&1 || {
        echo "[setup] installing Rust toolchain via rustup ..." >&2
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
        # shellcheck disable=SC1090
        . "$HOME/.cargo/env"
    }
    if [ -z "$SRC" ]; then
        SRC="$REFDIR/fastVEP"
        [ -d "$SRC/.git" ] || git clone https://github.com/Huang-lab/fastVEP.git "$SRC"
    fi
    echo "[setup] checking out fastVEP @ $FASTVEP_REF in $SRC ..." >&2
    ( cd "$SRC" && git fetch origin && git checkout "$FASTVEP_REF" )
    echo "[setup] building fastvep from $SRC ..." >&2
    ( cd "$SRC" && cargo install --path crates/fastvep-cli --force )
    command -v fastvep >/dev/null 2>&1 || {
        echo "ERROR: fastvep not on PATH after build — add ~/.cargo/bin to PATH" >&2; exit 3; }
    BUILT_FRESH=1
fi
echo "[setup] fastvep: $(fastvep --version 2>/dev/null || echo '(version unknown)')" >&2

# --- 1b. provenance file, so a run's exact fastvep build is on record ---
# Only written when this script actually built the binary against a known ref
# (see the WARN messages above for the two cases where provenance can't be
# established here). Append-only: each build adds a dated block rather than
# overwriting, so the file is a build history, not just the latest one.
if [ "$BUILT_FRESH" = 1 ]; then
    PROV="$REFDIR/fastvep.provenance.txt"
    RESOLVED_SHA="$(cd "$SRC" && git rev-parse HEAD)"
    {
        echo "=== fastvep build $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
        echo "fastvep_ref_requested: $FASTVEP_REF"
        echo "fastvep_commit_sha:    $RESOLVED_SHA"
        echo "fastvep_version:       $(fastvep --version 2>/dev/null || echo '(unknown)')"
        echo "source_repo:           https://github.com/Huang-lab/fastVEP"
        echo "source_path:           $SRC"
        echo "built_by:              ${USER:-unknown}"
        echo "built_on:              $(hostname 2>/dev/null || echo unknown)"
        echo
    } >> "$PROV"
    echo "[setup] provenance recorded: $PROV" >&2
fi

# --- 2. Ensembl GRCh38 GFF3 ---
GFF3="$REFDIR/Homo_sapiens.GRCh38.${REL}.gff3"
if [ -s "$GFF3" ]; then
    echo "[setup] GFF3 present: $GFF3" >&2
else
    url="https://ftp.ensembl.org/pub/release-${REL}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${REL}.gff3.gz"
    echo "[setup] downloading $url" >&2
    ( cd "$REFDIR" && curl -fSL -o "${GFF3}.gz" "$url" && gunzip -f "${GFF3}.gz" )
    echo "[setup] GFF3 ready: $GFF3" >&2
fi

# --- 3. FASTA index (reuse pipeline FASTA; only build .fai if missing) ---
if [ -n "$FASTA" ]; then
    if [ ! -s "${FASTA}.fai" ]; then
        if command -v samtools >/dev/null 2>&1; then
            echo "[setup] indexing FASTA (samtools faidx) ..." >&2
            samtools faidx "$FASTA"
        else
            echo "[setup] WARN: ${FASTA}.fai missing and samtools not found — index it before --hgvs runs" >&2
        fi
    fi
fi

# --- 4. optional transcript cache ---
if [ "$BUILD_CACHE" = 1 ]; then
    CACHE="$REFDIR/Homo_sapiens.GRCh38.${REL}.cache"
    if [ -s "$CACHE" ]; then
        echo "[setup] transcript cache present: $CACHE" >&2
    else
        echo "[setup] building transcript cache ..." >&2
        fastvep cache --gff3 "$GFF3" ${FASTA:+--fasta "$FASTA"} -o "$CACHE"
    fi
fi

# --- 5. print the env block to paste into runs ---
cat >&2 <<EOF

[setup] DONE. Use these in run_fastvep.sh / run_fastvep_batch.sh:

  export FASTVEP="$(command -v fastvep)"
  export GFF3="$GFF3"
${FASTA:+  export FASTA="$FASTA"}
${CACHE:+  export CACHE="$CACHE"}

  # smoke test (plain annotation, one chunk):
  acmg_fastvep/run_fastvep.sh -i <norm_qc/CHUNK.norm.vcf.gz> -o results-fastvep \\
      --gff3 "\$GFF3" ${FASTA:+--fasta "\$FASTA" --hgvs} ${CACHE:+--transcript-cache "\$CACHE"}
EOF
