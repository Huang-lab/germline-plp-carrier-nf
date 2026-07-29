# acmg_fastvep — standalone fastVEP annotation (and ACMG) on QC'd VCFs

Run **[fastVEP](https://github.com/Huang-lab/fastVEP)** (a Rust VEP reimplementation)
on the VCFs the Nextflow pipeline QC-normalizes — **outside** the pipeline, as plain
scripts. Two modes:

| Mode | Flag | Needs | Output |
|---|---|---|---|
| **Plain annotation** (default) | — | a GFF3 (+ FASTA for HGVS) | `<name>.fastvep.vcf` (CSQ: consequence, gene, HGVS, …) |
| **ACMG-AMP P/LP** | `--acmg` | GFF3 + FASTA + `--sa-dir` (supplementary DBs) | above **plus** `<name>.acmg_plp.tsv` (pipeline schema) |

**Why standalone:** it runs on the pipeline's already-published QC'd VCFs
(`results/<run>/norm_qc/<chunk>.norm.vcf.gz` — Ensembl contigs, PASS-only,
GT-masked), so fastVEP sees exactly the variants the ClinVar path does, without
touching `main.nf`. fastVEP is a native binary — **no Singularity/container or
`module load` needed** (unlike the Nextflow steps).

Plain annotation runs with **only a GFF3** — no gnomAD/ClinVar/REVEL downloads — so
you can prove fastVEP works on your data today, then layer ACMG on later.

## Files
- `setup_fastvep.sh` — one-time: build `fastvep` + fetch the Ensembl GRCh38 GFF3.
- `run_fastvep.sh` — annotate one VCF (plain, or `--acmg`). **Main entry point.**
- `run_fastvep_batch.sh` — run over a whole `norm_qc/` directory (local or one
  `bsub` per chunk); writes a chunk→status manifest.
- `run_fastvep_acmg.sh` — back-compat shim = `run_fastvep.sh --acmg`.
- `parse_fastvep_acmg.py` — CSQ `ACMG`/`ACMG_CRITERIA` → `acmg_plp.tsv` (used in
  `--acmg` mode). Reuses `plp_rules/csq.py`.

## 1. One-time setup
```bash
# reuse the SAME GRCh38 FASTA the pipeline uses (from params/msm.local.yaml)
acmg_fastvep/setup_fastvep.sh \
  --refdir /sc/arion/work/$USER/fastvep-refs \
  --fasta  /sc/arion/projects/rg_huangk06/variants_PLP_MSM/refs/vep_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# it prints the exact FASTVEP / GFF3 / FASTA exports to reuse below.
```
`--skip-build` if `fastvep` is already on PATH; `--build-cache` to pre-build the
transcript cache (otherwise auto-built next to the GFF3 on first run).

## 2. Plain annotation (the "just make it run" path)
One chunk:
```bash
export GFF3=/sc/arion/work/$USER/fastvep-refs/Homo_sapiens.GRCh38.115.gff3
export FASTA=/sc/arion/projects/rg_huangk06/variants_PLP_MSM/refs/vep_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa

acmg_fastvep/run_fastvep.sh \
  -i  results-chr1/norm_qc/<chunk>.norm.vcf.gz \
  -o  results-fastvep \
  --gff3 "$GFF3" --fasta "$FASTA" --hgvs
# -> results-fastvep/<chunk>.fastvep.vcf
```
All chunks in a run (local, sequential):
```bash
acmg_fastvep/run_fastvep_batch.sh \
  --in-dir results-chr1/norm_qc -o results-fastvep \
  --gff3 "$GFF3" --fasta "$FASTA" --hgvs
```
Or one LSF job per chunk (fastVEP is multithreaded):
```bash
acmg_fastvep/run_fastvep_batch.sh \
  --in-dir results-chr1/norm_qc -o results-fastvep \
  --gff3 "$GFF3" --fasta "$FASTA" --hgvs \
  --lsf -P "$MINERVA_ALLOCATION" --threads 4
# batch_manifest.tsv records chunk -> ok/FAILED/submitted (no silent drops).
```

## 3. ACMG mode (later — needs supplementary DBs)
Build the DBs once per fastVEP `docs/ACMG_SETUP.md` (gnomAD v4, ClinVar — pin the
same release as the ClinVar pipeline — REVEL, gene-level), each via
`fastvep sa-build --source <s> -i <download> -o <SA_DIR>/<s> --assembly GRCh38`
(`sa-build` is a **converter, not a downloader** — download the source first, or the
`.osa` will be empty). Then add `--acmg --sa-dir`:
```bash
acmg_fastvep/run_fastvep.sh \
  -i results-chr1/norm_qc/<chunk>.norm.vcf.gz -o results-fastvep-acmg \
  --gff3 "$GFF3" --fasta "$FASTA" --acmg --sa-dir /sc/arion/work/$USER/fastvep_sa
# -> <chunk>.fastvep.vcf  AND  <chunk>.acmg_plp.tsv
```
`acmg_plp.tsv` columns (identical to the pipeline's, so it is drop-in for the carrier
matrix later): `chr, pos, ref, alt, gene, acmg_label, acmg_criteria,
n_pathogenic_criteria, n_benign_criteria, is_acmg_PLP`. Per variant the parser keeps
the **most severe** call across transcripts (P > LP > VUS > LB > B); `is_acmg_PLP=1`
for P/LP.

## Notes
- **Contigs:** feed the pipeline's `*.norm.vcf.gz` (Ensembl `1,2,…`); the GFF3/FASTA
  (and any SA DBs) must be the same GRCh38/Ensembl space.
- The CSQ `ACMG`/`ACMG_CRITERIA` columns always appear in the header but are only
  *filled* when `--acmg` is passed.
- **Later integration:** because `acmg_plp.tsv` matches the pipeline schema, folding
  this into `main.nf` = a module that runs `run_fastvep.sh --acmg` and feeds the TSV
  into `CARRIER_GT`/`CARRIER_MATRIX` alongside the ClinVar table.

## Tests
```bash
python -m pytest acmg_fastvep/tests/ -q
```
- `test_parse_fastvep_acmg.py` — parser (severity collapse, criteria counts, P/LP
  flag) on a synthetic VCF; no fastVEP needed.
- `test_run_fastvep_smoke.py` — drives `run_fastvep.sh` on fastVEP's own
  `tests/test.vcf`; **auto-skips** when `fastvep` isn't on PATH (set `FASTVEP_SRC` if
  the checkout isn't `/workspace/fastvep`).
