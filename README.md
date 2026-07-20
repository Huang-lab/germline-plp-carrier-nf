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

## Quick start (local, no real data required)
```bash
# Unit tests
pip install pytest
python -m pytest tests/ -q

# DAG smoke test on synthetic fixtures
nextflow run . -profile test
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
