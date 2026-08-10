# germline-plp-carrier-nf

Cohort-agnostic Nextflow (DSL2) pipeline that turns WES pVCF chunks into
per-person carrier matrices under three P/LP definitions.

**Flow:** `pVCF chunks → norm+QC → VEP (ClinVar, AlphaMissense, LOFTEE, gnomAD)
→ P/LP sets → per-person carrier matrices`. ACMG-AMP classification is produced
by a **standalone [fastVEP](https://github.com/Huang-lab/fastVEP) tool
(`acmg_fastvep/`)** run on the QC'd VCFs, outside the Nextflow pipeline — see
[ACMG-AMP via fastVEP](#acmg-amp-via-fastvep).

## Design constraints
- **Code-only repo.** No participant data, no reference data, no results are
  ever tracked. A pre-commit hook enforces this — see `hooks/pre-commit`.
- **Genome-wide, all genes.** The cancer-gene panel is reporting/QC only and
  must never filter variants.
- **Gene symbol everywhere.** Every downstream table carries VEP `SYMBOL`.
- **ACMG-AMP via fastVEP.** ACMG P/LP is produced by the standalone
  `acmg_fastvep/` tool (Huang-lab/fastVEP), not by the Nextflow pipeline and not
  vendored — it is built on Minerva and run on the QC'd VCFs. (ANNOVAR/InterVar
  is no longer used.)

## Processing details

### Normalization & QC (`NORM_QC`)
For each input chunk, in order:
1. **Contig-naming reconciliation.** Input pVCFs are often UCSC-style
   (`chr1…chr22, chrX, chrY, chrM`), while the reference FASTA, VEP cache, and
   NCBI ClinVar are Ensembl-style (`1…22, X, Y, MT`). NORM_QC reads the
   chromosomes from the records and applies `bcftools annotate --rename-chrs`
   (`chrN → N`, `chrM/chrMT → MT`) so everything runs in one naming space.
   VCFs already in Ensembl style are passed through unchanged. **All outputs
   are therefore Ensembl-style contigs (`10`, not `chr10`).**
2. **Left-align + split multiallelics:** `bcftools norm -m -any -f <ref>`.
3. **Site filter:** `bcftools view -f PASS,.` — keeps only `PASS` sites
   (and filter-less `.`), dropping whatever the upstream caller flagged
   (e.g. `lowQual`). This respects any site-level QC already applied to the
   input.
4. **Per-genotype masking:** `bcftools +setGT … -i 'FMT/DP<${min_dp} | FMT/GQ<${min_gq}'`
   sets individual genotypes below the depth/quality thresholds to `./.`.
   This is finer-grained than a site-level average filter — it removes
   low-confidence *carrier calls* so a carrier is never called off a
   low-depth genotype. Defaults: `min_dp=10`, `min_gq=20`.

> **Note:** the allele-balance params (`het_ab_min/max`, `hom_ab_min`) are
> defined but **not yet applied** in `NORM_QC` (only DP and GQ are). Ask if you
> want het/hom AB masking wired into the `+setGT` expression.

### Classification
See the three P/LP frameworks below and `docs/data_dictionary.md`.

## Choosing classifiers
`params.classifiers` picks which in-pipeline P/LP definitions run — `clinvar`
and/or `am` (AlphaMissense). ACMG is **not** an in-pipeline classifier anymore;
it is produced by the standalone `acmg_fastvep/` tool (see below). Comma-separated
string or Groovy list. Examples:

```bash
nextflow run . -profile minerva -params-file params/msm.yaml --classifiers clinvar
nextflow run . -profile minerva -params-file params/msm.yaml \
    -params-file params/clinvar_only.yaml
```

`VALIDATE_CHUNK` only enforces CSQ subfields required by the selected
classifiers, so a ClinVar-only run does not need `am_pathogenicity` or gnomAD
fields in the annotated VCF.

## ACMG-AMP via fastVEP
ACMG-AMP P/LP classification is done with **[fastVEP](https://github.com/Huang-lab/fastVEP)**
(a Rust VEP reimplementation with a native `--acmg` flag: Richards 2015 +
ClinGen SVI), **not** ANNOVAR/InterVar. It lives in `acmg_fastvep/` as plain
scripts that run **outside** the Nextflow pipeline, on the QC'd VCFs the pipeline
publishes (`results/<run>/norm_qc/<chunk>.norm.vcf.gz`) — so ACMG is applied to
exactly the variants the ClinVar path sees, without destabilizing the pipeline.

Two modes:
- **Plain annotation** (default) — needs only an Ensembl GFF3 (+FASTA for HGVS);
  no supplementary databases. Fastest way to confirm fastVEP runs on your data.
- **ACMG P/LP** (`--acmg`) — adds `--sa-dir` supplementary DBs (gnomAD, ClinVar,
  REVEL, gene-level, built via `fastvep sa-build`) and emits `<chunk>.acmg_plp.tsv`
  in the same schema the carrier matrix consumes.

```bash
# one-time: build fastvep + fetch the GFF3 (reuses the pipeline's GRCh38 FASTA)
acmg_fastvep/setup_fastvep.sh --refdir /sc/arion/work/$USER/fastvep-refs --fasta <ref.fa>

# plain annotation over a whole run's QC output (local, or --lsf for one job/chunk)
acmg_fastvep/run_fastvep_batch.sh --in-dir results-chr1/norm_qc -o results-fastvep \
    --gff3 <...115.gff3> --fasta <ref.fa> --hgvs

# ACMG P/LP later, once supplementary DBs are built:
acmg_fastvep/run_fastvep.sh -i <chunk>.norm.vcf.gz -o results-fastvep-acmg \
    --gff3 <...115.gff3> --fasta <ref.fa> --acmg --sa-dir <sa_databases/>
```
Full details, setup, and tests: `acmg_fastvep/README.md`.

## Quick start (local, no real data required)
```bash
# Unit tests
pip install pytest
python -m pytest tests/ -q

# DAG smoke test on synthetic fixtures (Nextflow)
nextflow run . -profile test

# Same, without Nextflow (script emulates each process — used in authoring env)
tests/synthetic/run_dag_smoke.sh            # all three classifiers
tests/synthetic/run_clinvar_only_smoke.sh   # ClinVar only, chr21 fixture
```

## Real runs (Minerva)
See `docs/RUNNING_ON_MINERVA.md`.

## Outputs
- `results/carriers/carrier_matrix.tsv` (long; optional wide pivot)
- `results/variants/{clinvar_plp,acmg_plp,am_plp,qc_per_gene}.tsv`
- `results/manifest.json` — reference build, VEP + cache version, ClinVar
  release, AlphaMissense + calibration version, gnomAD version, pipeline git
  SHA, container digests.

Schemas: `docs/data_dictionary.md`, `assets/schemas/`.

## First-commit governance
```bash
ln -s ../../hooks/pre-commit .git/hooks/pre-commit
git ls-files | grep -E '\.(vcf|bcf|cram|bam)|SINAI_|SINAI-Million_' | grep -v '^tests/synthetic/' && echo BAD || echo OK
```
