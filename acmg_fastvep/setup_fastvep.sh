#!/usr/bin/env bash
# One-time setup of a pinned fastVEP (Huang-lab/fastVEP) install for the
# standalone ACMG path. Run on a Minerva LOGIN node (it downloads).
#
# Installs into one self-contained, versioned directory:
#
#   $FV_HOME/bin/fastvep          the v0.5.0 release binary (or a --from-source build)
#   $FV_HOME/src/                 fastVEP checked out at the v0.5.0 tag: the helper
#                                 scripts build_sa_dbs.sh / check_sa_db.sh call
#                                 (extract_gnomad_scores.py, check_acmg_stack.py,
#                                 clingen_gdv_to_oga.py, repeatmasker_to_bed.py)
#   $FV_HOME/ref/                 Ensembl 115 GFF3 + the fastVEP transcript cache
#   $FV_HOME/sa_src/  sa_db/      filled later by fetch_sa_sources.sh / build_sa_dbs.sh
#   $FV_HOME/fastvep.provenance   what was installed, from where, with which hashes
#   $FV_HOME/env.sh               `source` this in every job: FV_BIN, FV_CACHE, FV_SA ...
#
# Pinning. Everything is pinned to the v0.5.0 tag by commit SHA and to the
# release tarball by sha256. A mismatch on either is an error, never a warning:
# two runs months apart must not classify the same variant differently because
# a different binary happened to be found (Huang-lab/germline-plp-carrier-nf#3).
# The old reuse-whatever-is-on-PATH and --skip-build behaviours are gone for
# the same reason. Bump the three FASTVEP_* values below deliberately, in a
# commit that says why.
#
# The transcript cache needs ~16 GB of RAM to build, which is more than a login
# node should be asked for. By default this script prints the bsub command that
# builds it (--cache-only) instead of building it in place; pass --build-cache to
# build it here (e.g. when already inside an interactive compute session).
set -euo pipefail

FASTVEP_VERSION="v0.5.0"
FASTVEP_COMMIT="19a159f94fdb3d2aeac51812209387509e26814a"
FASTVEP_TARGET="x86_64-unknown-linux-gnu"
FASTVEP_TARBALL="fastVEP-${FASTVEP_VERSION}-${FASTVEP_TARGET}.tar.gz"
FASTVEP_TARBALL_SHA256="3b081b9e64711adbc0c2df73560df7d9a7240db05dc442ca811ddd8aa881ce82"
FASTVEP_REPO="https://github.com/Huang-lab/fastVEP"
FASTVEP_TARBALL_URL="$FASTVEP_REPO/releases/download/$FASTVEP_VERSION/$FASTVEP_TARBALL"
# The release workflow (.github/workflows/release.yml at the tag) packs
# "fastVEP-<tag>-<target>/fastvep"; extraction still searches if that moves.
FASTVEP_TARBALL_BIN="fastVEP-${FASTVEP_VERSION}-${FASTVEP_TARGET}/fastvep"

ENSEMBL_RELEASE="115"
GFF3_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz"

DEFAULT_FV_HOME="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/variants_PLP/data/fastvep-${FASTVEP_VERSION}"

usage() {
    cat >&2 <<USAGE
Usage: setup_fastvep.sh --fasta <GRCh38.fa> [options]

Required:
  --fasta <fa>             GRCh38 primary-assembly FASTA (Ensembl contig names, the
                           one the pipeline already uses). Not downloaded. A missing
                           .fai is built with samtools faidx.

Options:
  --fv-home <dir>          install root (default: \$FV_HOME, else
                           $DEFAULT_FV_HOME)
  --from-source            build the binary from the tag with cargo instead of using
                           the release tarball (cargo must be on PATH)
  --build-cache            build the transcript cache now (needs ~16 GB RAM)
  --cache-only             only build the transcript cache (binary, GFF3 and FASTA
                           must already be in place); what the printed bsub runs
  --allow-path-mismatch    do not fail when a different \`fastvep\` is first on PATH
  -h, --help

Pinned: fastVEP $FASTVEP_VERSION ($FASTVEP_COMMIT), Ensembl $ENSEMBL_RELEASE GFF3.
USAGE
    exit 2
}

FV_HOME="${FV_HOME:-$DEFAULT_FV_HOME}"
FASTA=""; FROM_SOURCE=0; BUILD_CACHE=0; CACHE_ONLY=0; ALLOW_PATH_MISMATCH=0
while [ $# -gt 0 ]; do
    case "$1" in
        --fasta)               FASTA="$2"; shift 2 ;;
        --fv-home)             FV_HOME="$2"; shift 2 ;;
        --from-source)         FROM_SOURCE=1; shift ;;
        --build-cache)         BUILD_CACHE=1; shift ;;
        --cache-only)          CACHE_ONLY=1; BUILD_CACHE=1; shift ;;
        --allow-path-mismatch) ALLOW_PATH_MISMATCH=1; shift ;;
        -h|--help)             usage ;;
        *) echo "Unknown arg: $1" >&2; usage ;;
    esac
done
[ -n "$FASTA" ] || { echo "ERROR: --fasta is required" >&2; usage; }
[ -s "$FASTA" ] || { echo "ERROR: FASTA missing or empty: $FASTA" >&2; exit 3; }
FASTA="$(cd "$(dirname "$FASTA")" && pwd)/$(basename "$FASTA")"
mkdir -p "$FV_HOME"
FV_HOME="$(cd "$FV_HOME" && pwd)"

BIN_DIR="$FV_HOME/bin"; SRC_DIR="$FV_HOME/src"; REF_DIR="$FV_HOME/ref"
DL_DIR="$FV_HOME/download"
FV_BIN="$BIN_DIR/fastvep"
GFF3="$REF_DIR/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3"
CACHE="$REF_DIR/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.fastvep.cache"
mkdir -p "$BIN_DIR" "$REF_DIR" "$DL_DIR" "$FV_HOME/sa_src" "$FV_HOME/sa_db"

sha256_of() {   # file -> hex digest (sha256sum on Linux, shasum elsewhere)
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}
fetch() {   # url dest — resumable, retried, fails on HTTP errors
    echo "[setup] downloading $1" >&2
    curl -fL -C - --retry 3 --retry-delay 10 -o "$2" "$1"
}

# --- 1. source tree at the tag (helper scripts; also the --from-source build) ---
if [ "$CACHE_ONLY" = 0 ]; then
    if [ ! -d "$SRC_DIR/.git" ]; then
        echo "[setup] cloning $FASTVEP_REPO @ $FASTVEP_VERSION -> $SRC_DIR" >&2
        git clone --quiet --branch "$FASTVEP_VERSION" "$FASTVEP_REPO.git" "$SRC_DIR"
    fi
fi
[ -d "$SRC_DIR/.git" ] || { echo "ERROR: no fastVEP source tree at $SRC_DIR (run without --cache-only first)" >&2; exit 3; }
SRC_COMMIT="$(git -C "$SRC_DIR" rev-parse HEAD)"
if [ "$SRC_COMMIT" != "$FASTVEP_COMMIT" ]; then
    echo "ERROR: $SRC_DIR is at $SRC_COMMIT, expected $FASTVEP_VERSION = $FASTVEP_COMMIT." >&2
    echo "       Remove it (or check out the tag) and re-run; the helper scripts must match the binary." >&2
    exit 4
fi
echo "[setup] source tree OK: $SRC_DIR @ $SRC_COMMIT" >&2

# --- 2. the binary ---
check_version() {   # bin -> 0 if it reports the pinned version
    local v
    v="$("$1" --version 2>/dev/null || true)"
    case "$v" in *"${FASTVEP_VERSION#v}"*) return 0 ;; *) return 1 ;; esac
}
if [ -x "$FV_BIN" ]; then
    check_version "$FV_BIN" || {
        echo "ERROR: $FV_BIN exists but reports '$("$FV_BIN" --version 2>&1 || true)', not $FASTVEP_VERSION." >&2
        echo "       Refusing to overwrite it; move it aside and re-run." >&2
        exit 4
    }
    echo "[setup] binary present: $FV_BIN ($("$FV_BIN" --version))" >&2
elif [ "$CACHE_ONLY" = 1 ]; then
    echo "ERROR: --cache-only but no binary at $FV_BIN" >&2; exit 3
elif [ "$FROM_SOURCE" = 1 ]; then
    command -v cargo >/dev/null 2>&1 || {
        echo "ERROR: --from-source needs cargo on PATH (rustup, or the Minerva rust module)" >&2; exit 3; }
    echo "[setup] building fastvep from $SRC_DIR (cargo build --release --locked) ..." >&2
    ( cd "$SRC_DIR" && cargo build --release --locked --bin fastvep )
    install -m 0755 "$SRC_DIR/target/release/fastvep" "$FV_BIN.partial"
    mv -f "$FV_BIN.partial" "$FV_BIN"
    INSTALL_METHOD="cargo build --release --locked (from $SRC_DIR)"
else
    TB="$DL_DIR/$FASTVEP_TARBALL"
    if [ ! -s "$TB" ] || [ "$(sha256_of "$TB")" != "$FASTVEP_TARBALL_SHA256" ]; then
        fetch "$FASTVEP_TARBALL_URL" "$TB"
    fi
    got="$(sha256_of "$TB")"
    if [ "$got" != "$FASTVEP_TARBALL_SHA256" ]; then
        echo "ERROR: $TB sha256 $got != pinned $FASTVEP_TARBALL_SHA256. Not installing." >&2
        echo "       Delete it and re-run; if it repeats, the release asset changed upstream." >&2
        exit 4
    fi
    echo "[setup] tarball sha256 OK ($got)" >&2
    XT="$(mktemp -d "$DL_DIR/extract.XXXXXX")"
    tar -xzf "$TB" -C "$XT"
    cand="$XT/$FASTVEP_TARBALL_BIN"
    if [ ! -f "$cand" ]; then
        # Layout moved: accept exactly one regular file named fastvep anywhere in it.
        n="$(find "$XT" -type f -name fastvep | wc -l | tr -d ' ')"
        [ "$n" = 1 ] || { echo "ERROR: expected one 'fastvep' in $TB, found $n" >&2; rm -rf "$XT"; exit 4; }
        cand="$(find "$XT" -type f -name fastvep)"
        echo "[setup] WARN: binary not at $FASTVEP_TARBALL_BIN; using $cand" >&2
    fi
    install -m 0755 "$cand" "$FV_BIN.partial"
    mv -f "$FV_BIN.partial" "$FV_BIN"
    rm -rf "$XT"
    INSTALL_METHOD="release tarball $FASTVEP_TARBALL_URL"
fi
check_version "$FV_BIN" || {
    echo "ERROR: $FV_BIN reports '$("$FV_BIN" --version 2>&1 || true)', not $FASTVEP_VERSION" >&2
    echo "       (a glibc error here means this host is older than the release supports; try --from-source)" >&2
    exit 4
}
FV_BIN_SHA256="$(sha256_of "$FV_BIN")"
FV_VERSION_STR="$("$FV_BIN" --version)"
echo "[setup] binary: $FV_BIN  $FV_VERSION_STR  sha256=$FV_BIN_SHA256" >&2

# A different fastvep first on PATH (the old 0.3.0 install is) is how a job
# that forgets FV_BIN would silently run the wrong version. Every script here
# calls "$FV_BIN" explicitly, but say so and stop unless told it is understood.
if path_bin="$(command -v fastvep 2>/dev/null)" && [ -n "$path_bin" ]; then
    path_sha="$(sha256_of "$path_bin")"
    if [ "$path_sha" != "$FV_BIN_SHA256" ]; then
        echo "[setup] CONFLICT: \`fastvep\` on PATH is $path_bin" >&2
        echo "          ($("$path_bin" --version 2>/dev/null || echo 'version unknown'), sha256=$path_sha)" >&2
        echo "        and is NOT the pinned $FV_BIN (sha256=$FV_BIN_SHA256)." >&2
        echo "        The acmg_fastvep scripts call \$FV_BIN from env.sh explicitly and never" >&2
        echo "        the PATH one, but anything typed by hand as plain \`fastvep\` gets the old" >&2
        echo "        binary. Remove it from PATH, or pass --allow-path-mismatch." >&2
        [ "$ALLOW_PATH_MISMATCH" = 1 ] || exit 5
        echo "[setup] --allow-path-mismatch: continuing" >&2
    fi
fi

# --- 3. Ensembl GFF3 ---
if [ -s "$GFF3" ]; then
    echo "[setup] GFF3 present: $GFF3" >&2
elif [ "$CACHE_ONLY" = 1 ]; then
    echo "ERROR: --cache-only but no GFF3 at $GFF3" >&2; exit 3
else
    fetch "$GFF3_URL" "$DL_DIR/$(basename "$GFF3").gz"
    gzip -t "$DL_DIR/$(basename "$GFF3").gz"
    gzip -dc "$DL_DIR/$(basename "$GFF3").gz" > "$GFF3.partial"
    mv -f "$GFF3.partial" "$GFF3"
    echo "[setup] GFF3 ready: $GFF3" >&2
fi

# --- 4. FASTA index ---
if [ ! -s "${FASTA}.fai" ]; then
    command -v samtools >/dev/null 2>&1 || {
        echo "ERROR: ${FASTA}.fai missing and samtools not on PATH (module load samtools)" >&2; exit 3; }
    echo "[setup] samtools faidx $FASTA" >&2
    samtools faidx "$FASTA"
fi

# --- 5. transcript cache (built by THIS binary; a cache from 0.3.0 is not reused) ---
CACHE_LOG="$REF_DIR/$(basename "$CACHE").build.log"
if [ -s "$CACHE" ]; then
    echo "[setup] transcript cache present: $CACHE" >&2
elif [ "$BUILD_CACHE" = 1 ]; then
    echo "[setup] building transcript cache (~16 GB RAM) -> $CACHE" >&2
    rm -f "$CACHE.partial"
    "$FV_BIN" cache --gff3 "$GFF3" --fasta "$FASTA" -o "$CACHE.partial" --no-progress 2> "$CACHE_LOG" || {
        cat "$CACHE_LOG" >&2; rm -f "$CACHE.partial"; echo "ERROR: fastvep cache failed" >&2; exit 6; }
    [ -s "$CACHE.partial" ] || { echo "ERROR: fastvep cache wrote nothing" >&2; exit 6; }
    mv -f "$CACHE.partial" "$CACHE"
    echo "[setup] transcript cache ready: $CACHE" >&2
fi

# --- 6. provenance + env.sh (rewritten on every run, so they describe what is there now) ---
glibc="$( (getconf GNU_LIBC_VERSION 2>/dev/null || ldd --version 2>/dev/null | head -1) || echo unknown)"
cache_sha=""; [ -s "$CACHE" ] && cache_sha="$(sha256_of "$CACHE")"
PROV="$FV_HOME/fastvep.provenance"
{
    echo "# fastVEP install provenance, written by setup_fastvep.sh"
    echo "fastvep_tag              $FASTVEP_VERSION"
    echo "fastvep_commit           $FASTVEP_COMMIT"
    echo "source_tree_commit       $SRC_COMMIT"
    echo "install_method           ${INSTALL_METHOD:-existing $FV_BIN (installed by an earlier run)}"
    echo "release_tarball          $FASTVEP_TARBALL"
    echo "release_tarball_sha256   $FASTVEP_TARBALL_SHA256"
    echo "binary                   $FV_BIN"
    echo "binary_sha256            $FV_BIN_SHA256"
    echo "binary_version           $FV_VERSION_STR"
    echo "gff3                     $GFF3"
    echo "gff3_url                 $GFF3_URL"
    echo "fasta                    $FASTA"
    echo "transcript_cache         $CACHE"
    echo "transcript_cache_sha256  ${cache_sha:-NOT_BUILT}"
    echo "glibc                    $glibc"
    echo "date_utc                 $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "host                     $(hostname 2>/dev/null || echo unknown)"
    echo "user                     ${USER:-unknown}"
} > "$PROV"

ENV_SH="$FV_HOME/env.sh"
cat > "$ENV_SH" <<EOF
# Written by acmg_fastvep/setup_fastvep.sh ($(date -u +%Y-%m-%dT%H:%M:%SZ)). Source me:
#   source $ENV_SH
export FV_HOME="$FV_HOME"
export FV_BIN="$FV_BIN"
export FV_BIN_SHA256="$FV_BIN_SHA256"
export FV_GFF3="$GFF3"
export FV_FASTA="$FASTA"
export FV_CACHE="$CACHE"
export FV_SA="$FV_HOME/sa_db"
export FV_SRC="$SRC_DIR"
EOF
echo "[setup] wrote $PROV and $ENV_SH" >&2

if [ ! -s "$CACHE" ]; then
    HERE="$(cd "$(dirname "$0")" && pwd)"
    cat >&2 <<EOF

[setup] Transcript cache NOT built yet (needs ~16 GB RAM; not done on a login node).
        Build it as a job, then re-run nothing else:

  mkdir -p logs
  bsub -P "\$MINERVA_ALLOCATION" -q premium -n 1 -R "rusage[mem=20000]" -W 2:00 \\
       -J fv050-cache -oo logs/fv050-cache.%J.out -eo logs/fv050-cache.%J.err \\
       "$HERE/setup_fastvep.sh --fv-home $FV_HOME --fasta $FASTA --cache-only --allow-path-mismatch"
EOF
fi

cat >&2 <<EOF

[setup] DONE. Next:
  source $ENV_SH
  acmg_fastvep/fetch_sa_sources.sh            # login node: downloads into \$FV_HOME/sa_src
  bsub < acmg_fastvep/build_sa_dbs.lsf         # see acmg_fastvep/README.md
EOF
