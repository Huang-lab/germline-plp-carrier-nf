# acmg_fastvep — standalone ACMG-AMP via fastVEP

ACMG-AMP P/LP classification using **[fastVEP](https://github.com/Huang-lab/fastVEP)**
(`--acmg`: full Richards 2015 + ClinGen SVI, 28 criteria) instead of
ANNOVAR/InterVar.

**This is intentionally OUTSIDE the Nextflow pipeline.** It runs after the
pipeline's QC step, on the QC'd VCFs the pipeline already publishes, so ACMG is
applied to exactly the same variants as the ClinVar calls — without touching the
working ClinVar pipeline. Once validated, folding it into `main.nf` (replacing
the ANNOVAR/InterVar module) is a small, low-risk step.

## What it does
```
QC'd VCF (results/<run>/norm_qc/*.norm.vcf)   # Ensembl contigs, PASS-only, GT-masked
      │
      ▼  fastvep annotate --acmg --output-format vcf
<name>.fastvep.vcf                            # CSQ incl. ACMG + ACMG_CRITERIA
      │
      ▼  parse_fastvep_acmg.py
<name>.acmg_plp.tsv                           # pipeline acmg_plp.tsv schema
```
Output columns (identical to the pipeline's `acmg_plp.tsv`, so it's drop-in for
the carrier matrix later):
`chr, pos, ref, alt, gene, acmg_label, acmg_criteria, n_pathogenic_criteria, n_benign_criteria, is_acmg_PLP`

Per-variant collapse: fastVEP emits per-transcript ACMG; the parser keeps the
**most severe** call across transcripts (P > LP > VUS > LB > B). `is_acmg_PLP` = 1
for P or LP.

## One-time setup on Minerva

### 1. Install fastVEP
```bash
# Rust toolchain, then build
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh && source "$HOME/.cargo/env"
git clone https://github.com/Huang-lab/fastVEP.git && cd fastVEP
cargo install --path crates/fastvep-cli     # installs `fastvep` to ~/.cargo/bin
fastvep --version
```
(Or the conda recipe under `fastVEP/conda/recipe/` — see fastVEP README.)

### 2. Gene models + reference
```bash
# GFF3 (Ensembl) and the SAME GRCh38 FASTA the pipeline uses (bgzipped+faidx ok)
#   Homo_sapiens.GRCh38.<rel>.gff3
#   Homo_sapiens.GRCh38.dna.primary_assembly.fa(.bgz) + .fai
```

### 3. Supplementary databases (`$SA_DIR`, the big step)
Follow **fastVEP `docs/ACMG_SETUP.md`**. Minimum viable ACMG (~90% of criteria):
- **gnomAD v4** (BA1/BS1/BS2/PM2) — largest source
- **ClinVar** (PP5/BP6 optional; PM2 frequency backstop) — pin the **same
  release date** you used for the ClinVar pipeline for consistency
- **REVEL** (PP3/BP4 missense)
- gene-level `.oga` (PVS1/PS1/PM1/PM5 — constraint, ClinGen GDV)

Each is `fastvep sa-build --source <s> -i <downloaded> -o $SA_DIR/<s> --assembly GRCh38`.
Verify `.osa` sizes are non-trivial (a few-KB `.osa` = empty build).

## Run it (per chunk, or loop over the pipeline's QC output)
```bash
export SA_DIR=/sc/arion/work/$USER/fastvep_sa
GFF3=/sc/arion/work/$USER/germline-plp-refs/Homo_sapiens.GRCh38.115.gff3
FASTA=/sc/arion/projects/rg_huangk06/variants_PLP_MSM/refs/vep_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz

acmg_fastvep/run_fastvep_acmg.sh \
  -i  /sc/arion/projects/rg_huangk06/variants_PLP_MSM/results-chr1/norm_qc/<chunk>.norm.vcf \
  -o  /sc/arion/projects/rg_huangk06/variants_PLP_MSM/results-chr1-acmg \
  --gff3 "$GFF3" --fasta "$FASTA" --sa-dir "$SA_DIR"
```
Run under LSF for many chunks (loop, or one `bsub` per chunk) — fastVEP is
multi-threaded, so give it several cores.

## Notes
- **Contigs:** feed the pipeline's `*.norm.vcf` (already stripped to Ensembl
  `1,2,…`). Ensure your GFF3/FASTA/SA are the same GRCh38/Ensembl space.
- **Thresholds:** override any ACMG threshold via `--acmg-config <config.toml>`
  (see fastVEP `docs/ACMG.md`).
- **Later integration:** because the output matches `acmg_plp.tsv`, wiring this
  into the pipeline = a new module that runs `run_fastvep_acmg.sh` and feeds the
  TSV into `CARRIER_GT`/`CARRIER_MATRIX` alongside the ClinVar table.

## Tests
```bash
python -m pytest acmg_fastvep/tests/ -q
```
Validates the CSQ→acmg_plp parser (severity collapse, criteria counts,
P/LP flag) on a synthetic fastVEP VCF — no fastVEP install needed.
