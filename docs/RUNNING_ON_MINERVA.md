# Running germline-plp-carrier-nf on Minerva

## 1. Load modules & set up Nextflow
Nextflow is installed via conda, not the site `nextflow` module (the site
module version drifts).

```bash
ml java anaconda3 singularity-ce
conda create -y -n nextflow -c bioconda nextflow=23.10.*
conda activate nextflow
```

## 2. Proxy (required for Singularity pulls, curl, etc.)
```bash
export http_proxy=http://172.28.7.1:3128
export https_proxy=http://172.28.7.1:3128
export all_proxy=http://172.28.7.1:3128
export no_proxy="*.chimera.hpc.mssm.edu,.mssm.edu,localhost,127.0.0.1"
```

## 3. LSF allocation
```bash
export MINERVA_ALLOCATION=acc_<your_project>
```

## 4. Obtain resources (run scripts here, not in the authoring env)
```bash
export RESOURCES_DIR=/sc/arion/work/$USER/germline-plp-refs
export SHARED_REFS=/sc/arion/projects/YOUR_PROJECT/refs   # optional
export CLINVAR_RELEASE=20260401                            # pin one date
export VEP_CACHE_VERSION=113
export GNOMAD_VERSION=v4.1

setup/fetch_references.sh          # reuses SHARED_REFS/* when present

export ANNOVAR_TARBALL=/path/to/your/annovar.latest.tar.gz
setup/setup_annovar_intervar.sh    # SAME CLINVAR_RELEASE — critical
```

Then supply the AlphaMissense gene-specific calibration TSV (Chen/Pejaver
2026) yourself and set `params.am_calibration_tsv` to its path in
`params/<cohort>.yaml`.

## 5. Build the two Singularity images
See `containers/README.md`. Do **not** push either image to a public registry.

## 6. Verify the gnomAD popmax field name
Confirm the popmax AF field name in the real gnomAD v4 VCF header on Minerva:

```bash
zgrep -m1 '^##INFO=<ID=AF_grpmax' $RESOURCES_DIR/gnomad/gnomad.exomes.v4.1.sites.vcf.bgz
```

Adjust `params.gnomad_popmax_field` if the field is named differently in the
release you use.

## 7. Submit the Nextflow driver as an LSF job
Do NOT run the driver on a login node. Use `run_minerva.lsf`:

```bash
bsub -P $MINERVA_ALLOCATION < run_minerva.lsf
```

The driver fans each process out as its own child bsub. First validate a
single chunk before scaling:

```bash
# Single-chunk validation
nextflow run . -profile minerva -params-file params/msm.yaml \
    --input_vcfs '/sc/arion/projects/CHANGEME/msm/wes/pvcf_chunks/chunk_001.vcf.gz'
```

## 8. Output check
```bash
head results/carriers/carrier_matrix.tsv
jq . results/manifest.json
```
