# Running germline-plp-carrier-nf on Minerva

## 1. Load modules & set up Nextflow
Nextflow is installed via conda, not the site `nextflow` module (the site
module version drifts).

```bash
ml java anaconda3 singularity-ce
conda create -y -n nextflow -c bioconda nextflow=23.10.*
conda activate nextflow
```

## 2. Proxy (usually NOT needed)
Minerva login/compute nodes have direct internet access. **Do not set an HTTP
proxy** unless a login shows you're on a subnet that requires one — the
nf-core `mssm` config also sets no proxy.

If you've already exported bad values, unset them:
```bash
unset http_proxy https_proxy all_proxy no_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY NO_PROXY
curl -sSI https://repo.anaconda.com/ | head -1   # sanity: expect HTTP/2 200
```
If direct access fails on your node, ask HPC support (hpchelp@hpc.mssm.edu) for
the current proxy for your subnet, then export those values before conda /
`singularity build` / `setup/fetch_references.sh`.

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

## 5. Containers (biocontainers, no build required)
By default `conf/minerva.config` points each process at a **pre-built public
Docker image**, and Singularity pulls it on first use into
`$SINGULARITY_CACHEDIR`. No `singularity build` step is needed — this
sidesteps Minerva's no-fakeroot restriction.

Defaults:
- `withLabel: vep` → `docker://ensemblorg/ensembl-vep:release_113.0`
- `withLabel: bcf` → `docker://staphb/bcftools:1.20`
- `withLabel: py`  → `docker://python:3.11-slim`
- `withLabel: annovar` → user-supplied (only needed for ACMG classifier)

First pipeline run will pull these once (~10 min total, one-time cost).

**Optional:** If you'd rather build a custom image, `containers/annotate.def`
still exists; build it with `singularity build --fakeroot` (only works if
fakeroot is enabled for your account) or `--remote` (Sylabs cloud), then set
`params.container_annotate` to the .sif path — it overrides the biocontainer
default for `vep/bcf/py` labels.

## 6. Verify the gnomAD popmax field name
Confirm the popmax AF field name in the real gnomAD v4 VCF header on Minerva:

```bash
zgrep -m1 '^##INFO=<ID=AF_grpmax' $RESOURCES_DIR/gnomad/gnomad.exomes.v4.1.sites.vcf.bgz
```

Adjust `params.gnomad_popmax_field` if the field is named differently in the
release you use.

## 6a. (Optional) Run a single classifier
To validate the ClinVar-only path first (fastest; needs no ANNOVAR or
AlphaMissense data):

```bash
nextflow run . -profile minerva -params-file params/msm.yaml \
    --classifiers clinvar \
    --input_vcfs '/sc/arion/projects/CHANGEME/msm/wes/pvcf_chunks/chr21.*.vcf.gz'
```

Any subset of `clinvar,acmg,am` is valid.

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
