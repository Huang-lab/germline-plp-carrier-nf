# germline-plp-carrier-nf

Cohort-agnostic Nextflow (DSL2) pipeline that turns WES pVCF chunks into
per-person carrier matrices under three P/LP definitions.

**Flow:** `pVCF chunks → norm+QC → VEP (ClinVar, AlphaMissense, LOFTEE, gnomAD)
→ ANNOVAR+InterVar (ACMG) → three P/LP sets → per-person carrier matrices`.

## Design constraints
- **Code-only repo.** No participant data, no reference data, no results are
  ever tracked. A pre-commit hook enforces this — see `hooks/pre-commit`.
- **Genome-wide, all genes.** The cancer-gene panel is reporting/QC only and
  must never filter variants.
- **Gene symbol everywhere.** Every downstream table carries VEP `SYMBOL`.
- **ANNOVAR is user-supplied.** Not vendored. The container builds from your
  registration-gated tarball on Minerva and is never pushed publicly.

## Choosing classifiers
`params.classifiers` picks which P/LP definitions run. Any subset of
`clinvar`, `acmg`, `am`. Comma-separated string or Groovy list. Default: all
three. Examples:

```bash
nextflow run . -profile minerva -params-file params/msm.yaml --classifiers clinvar
nextflow run . -profile minerva -params-file params/msm.yaml \
    -params-file params/clinvar_only.yaml
```

`VALIDATE_CHUNK` only enforces CSQ subfields required by the selected
classifiers, so a ClinVar-only run does not need `am_pathogenicity` or gnomAD
fields in the annotated VCF.

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
