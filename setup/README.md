# setup/ — resource acquisition (Minerva only)

These scripts are **authored** in this repo but **run on Minerva**. Do not
execute them in the authoring environment — they download hundreds of GB and
must land on Minerva scratch/work.

## `fetch_references.sh` (automated)
Downloads / symlinks:
- ClinVar VCF+tbi (**pinned** by `$CLINVAR_RELEASE`)
- VEP cache + FASTA (`$VEP_CACHE_VERSION`)
- AlphaMissense TSV+tbi + Ensembl VEP AlphaMissense plugin
- LOFTEE plugin data (GRCh38: ancestral FASTA, conservation, GERP)
- gnomAD v4 exome sites (huge — **reuses `$SHARED_REFS/gnomad/…` if present**)

It writes `$RESOURCES_DIR/versions.txt` for the run manifest.

### Required env
```bash
export RESOURCES_DIR=/sc/arion/work/$USER/germline-plp-refs
export CLINVAR_RELEASE=20260401   # pick a dated release; MUST match ANNOVAR humandb below
export VEP_CACHE_VERSION=113
export GNOMAD_VERSION=v4.1
export SHARED_REFS=/sc/arion/projects/YOUR_PROJECT/refs   # optional; reuses existing lab paths
setup/fetch_references.sh
```

## `setup_annovar_intervar.sh` (license-gated)
ANNOVAR is registration-gated. You must obtain the tarball yourself and pass
its path via `$ANNOVAR_TARBALL`. The script fails clearly if the tarball is
missing — it **never fetches ANNOVAR**.

InterVar is freely available (MIT) and cloned by the script. The humandb DBs
InterVar needs are then downloaded via `annotate_variation.pl -downdb -webfrom annovar`.

### ClinVar-date pinning (critical)
The `clinvar_<date>` humandb version installed here MUST match the ClinVar VCF
release fetched by `fetch_references.sh` — pass the same `$CLINVAR_RELEASE`.
Otherwise the VEP `--custom` ClinVar track (used by ClinVar classifier) and the
ANNOVAR/InterVar ACMG track will disagree.

```bash
export ANNOVAR_TARBALL=/path/to/annovar.latest.tar.gz
export RESOURCES_DIR=/sc/arion/work/$USER/germline-plp-refs
export CLINVAR_RELEASE=20260401   # SAME as above
setup/setup_annovar_intervar.sh
```

### User-supplied (not downloaded)
- **AlphaMissense gene-specific calibration table** (Chen/Pejaver 2026): supply
  the TSV yourself and set `params.am_calibration_tsv` to its path.
- **ANNOVAR tarball**: as above.

## After setup
Edit `params/msm.yaml` or `params/biome.yaml` to point at the paths under
`$RESOURCES_DIR/` (and your two `.sif` images built per `containers/README.md`).
