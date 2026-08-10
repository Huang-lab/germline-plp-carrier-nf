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
export CLINVAR_RELEASE=20260401   # pick a dated release (recorded in the manifest)
export VEP_CACHE_VERSION=113
export GNOMAD_VERSION=v4.1
export SHARED_REFS=/sc/arion/projects/YOUR_PROJECT/refs   # optional; reuses existing lab paths
setup/fetch_references.sh
```

## ACMG via fastVEP (`--classifiers acmg`)
ACMG-AMP is produced by **fastVEP** (Huang-lab/fastVEP), not ANNOVAR/InterVar.
Setup is separate and documented in `acmg_fastvep/README.md`:
```bash
acmg_fastvep/setup_fastvep.sh --refdir /sc/arion/work/$USER/germline-plp-refs/fastvep \
    --fasta <the same GRCh38 FASTA used above>
```
Then build the supplementary databases (gnomAD/ClinVar/REVEL/gene-level) per
fastVEP `docs/ACMG_SETUP.md`, and set `params.fastvep_gff3` + `params.fastvep_sa_dir`
(and optionally `params.fastvep_bin` / `params.container_fastvep`). For
consistency, build the ClinVar SA DB from the SAME `$CLINVAR_RELEASE` used above.

### User-supplied (not downloaded)
- **AlphaMissense gene-specific calibration table** (Chen/Pejaver 2026): supply
  the TSV yourself and set `params.am_calibration_tsv` to its path.

## After setup
Edit `params/msm.yaml` or `params/biome.yaml` to point at the paths under
`$RESOURCES_DIR/` (and your `.sif` images pulled per `docs/RUNNING_ON_MINERVA.md`).
