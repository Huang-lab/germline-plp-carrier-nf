# Running germline-plp-carrier-nf on Minerva

## 1. Load modules & set up Nextflow
Nextflow is installed via conda, not the site `nextflow` module (the site
module version drifts).

```bash
ml java anaconda3 singularity-ce
conda create -y -n nextflow -c bioconda nextflow=23.10.*
conda activate nextflow
```

## 2. Proxy (REQUIRED)
Per the nf-core `mssm` doc, Minerva requires a proxy for all outbound network
(conda, container pulls, reference downloads). Export before any download step
and inside every bsub:
```bash
export http_proxy=http://172.28.7.1:3128
export https_proxy=http://172.28.7.1:3128
export all_proxy=http://172.28.7.1:3128
export no_proxy=localhost,*.chimera.hpc.mssm.edu,172.28.0.0/16
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
```

**ACMG (`--classifiers acmg`) uses fastVEP, not ANNOVAR/InterVar.** Set it up
separately (native binary + supplementary DBs) per `acmg_fastvep/README.md`
(`acmg_fastvep/setup_fastvep.sh`), then set `params.fastvep_gff3` +
`params.fastvep_sa_dir`. Not needed for ClinVar-only or AlphaMissense-only runs.

Then supply the AlphaMissense gene-specific calibration TSV (Chen/Pejaver
2026) yourself and set `params.am_calibration_tsv` to its path in
`params/<cohort>.yaml`.

## 5. Containers — pull public images ON A COMPUTE NODE
Minerva does not permit `singularity build` (no fakeroot). Instead, pull three
pre-built public images. **Login nodes cap per-user threads and crash during
OCI→SIF conversion** (`FATAL ERROR: Failed to create thread`), so pull from a
compute node via an interactive job:

```bash
bsub -P $MINERVA_ALLOCATION -q premium -n 4 -W 2:00 \
     -R "rusage[mem=16000] span[hosts=1]" -Is bash

# inside the interactive shell:
ml singularity-ce
export SINGULARITY_CACHEDIR=/sc/arion/work/$USER/singularity_cache
export SINGULARITY_TMPDIR=/sc/arion/work/$USER/singularity_tmp
mkdir -p $SINGULARITY_CACHEDIR $SINGULARITY_TMPDIR
export http_proxy=http://172.28.7.1:3128 https_proxy=http://172.28.7.1:3128 all_proxy=http://172.28.7.1:3128

singularity pull $SINGULARITY_CACHEDIR/ensembl-vep_release_113.0.sif docker://ensemblorg/ensembl-vep:release_113.0
singularity pull $SINGULARITY_CACHEDIR/bcftools_1.20.sif             docker://staphb/bcftools:1.20
singularity pull $SINGULARITY_CACHEDIR/python_3.11-slim.sif          docker://python:3.11-slim
ls -lh $SINGULARITY_CACHEDIR/*.sif
exit
```

Then point `params/<cohort>.yaml` at the absolute `.sif` paths:
```yaml
container_vep:      '/sc/arion/work/<USER>/singularity_cache/ensembl-vep_release_113.0.sif'
container_bcftools: '/sc/arion/work/<USER>/singularity_cache/bcftools_1.20.sif'
container_python:   '/sc/arion/work/<USER>/singularity_cache/python_3.11-slim.sif'
```
Referencing local `.sif` files (not `docker://` URLs) means Nextflow reuses the
images directly — no re-pull, no cache-name mismatch, works offline.

Process→image map: NORM_QC→bcftools, VEP_ANNOTATE→vep, everything else→python.
ACMG_FASTVEP runs the native `fastvep` binary (no container) unless you set
`params.container_fastvep`.

## 6. Verify the gnomAD popmax field name
Confirm the popmax AF field name in the real gnomAD v4 VCF header on Minerva:

```bash
zgrep -m1 '^##INFO=<ID=AF_grpmax' $RESOURCES_DIR/gnomad/gnomad.exomes.v4.1.sites.vcf.bgz
```

Adjust `params.gnomad_popmax_field` if the field is named differently in the
release you use.

## 6a. (Optional) Run a single classifier
To validate the ClinVar-only path first (fastest; needs no fastVEP or
AlphaMissense data):

```bash
nextflow run . -profile minerva -params-file params/msm.yaml \
    --classifiers clinvar \
    --input_vcfs '/sc/arion/projects/CHANGEME/msm/wes/pvcf_chunks/chr21.*.vcf.gz'
```

Any subset of `clinvar,acmg,am` is valid.

## 7. Submit the Nextflow driver as an LSF job
Do NOT run the driver on a login node. The driver's bsub MUST export the proxy
(so LSF child submission works) and `NXF_OFFLINE=true` (so Nextflow does not try
to reach nextflow.io on startup — that call fails on compute nodes and stalls
the run). Use `run_minerva.lsf`, or submit inline:

```bash
SMALLEST=$(ls -Sr $INPUT_DIR/*.vcf.gz | head -1)   # pilot: smallest chunk
bsub -P $MINERVA_ALLOCATION -q premium -n 2 -W 24:00 \
     -R "rusage[mem=8000] span[hosts=1]" \
     -J plp-pilot -oo logs/pilot.%J.out -eo logs/pilot.%J.err <<EOF
ml java anaconda3 singularity-ce
source activate nextflow
export MINERVA_ALLOCATION=$MINERVA_ALLOCATION
export SINGULARITY_CACHEDIR=/sc/arion/work/\$USER/singularity_cache
export SINGULARITY_TMPDIR=/sc/arion/work/\$USER/singularity_tmp
export http_proxy=http://172.28.7.1:3128 https_proxy=http://172.28.7.1:3128 all_proxy=http://172.28.7.1:3128
export no_proxy=localhost,*.chimera.hpc.mssm.edu,172.28.0.0/16
export NXF_OFFLINE=true
export NXF_DISABLE_CHECK_LATEST=true
cd $ROOT/germline-plp-carrier-nf
nextflow run . -profile minerva -params-file params/msm.yaml \\
    --classifiers clinvar --input_vcfs "\$SMALLEST" --outdir "$ROOT/results-pilot"
EOF
```

After the pilot succeeds, run all chunks: drop `--input_vcfs`/`--outdir` (they
come from the params file) and add `-resume`.

## 8. Output check
```bash
head results/carriers/carrier_matrix.tsv
jq . results/manifest.json
```
