# Data dictionary

All outputs are **de-identified per site policy**: a `person_id` column
contains only the site-native ID. This pipeline does not link to phenotypes.

## `results/carriers/carrier_matrix.tsv` (long)

| column          | type   | description |
|-----------------|--------|-------------|
| `chr`           | str    | Chromosome (e.g. `chr17`). |
| `pos`           | int    | 1-based position (post-`bcftools norm`, left-aligned, biallelic). |
| `ref`, `alt`    | str    | REF and ALT alleles. |
| `gene`          | str    | VEP `SYMBOL`. |
| `person_id`     | str    | Sample ID as it appears in the input VCF; filtered by `params.keep_samples` when supplied. |
| `GT`            | str    | Raw genotype as called (e.g. `0/1`, `1/1`, `1`). |
| `zygosity`      | str    | Derived from `GT`: `het` / `hom_alt` / `hemizygous` (haploid call, or single allele on chrX/Y). |
| `is_clinvar_PLP`| 0/1    | 1 iff variant is P/LP per ClinVar with review-status ≥ `params.clinvar_min_stars`. |
| `is_acmg_PLP`   | 0/1    | 1 iff fastVEP's ACMG-AMP call is Pathogenic or Likely pathogenic. |
| `is_AM_PLP`     | 0/1    | 1 iff AlphaMissense score meets the gene-specific threshold at `params.am_min_strength`. |

Rows exist only for carriers (at least one alt allele in the person's genotype)
whose variant is qualifying under at least one of the three definitions.
Contigs are Ensembl-style (e.g. `17`, not `chr17`) — NORM_QC strips the `chr`
prefix so the data matches the Ensembl reference / VEP cache / ClinVar.

## `results/carriers/carrier_matrix.wide.tsv` (optional pivot)
Rows: `variant_key = chr:pos:ref:alt`. Columns: samples. Cell = 1 iff carrier
of a variant that is P/LP under any framework, else 0.

## `results/variants/clinvar_plp.tsv`
| column | description |
|--------|-------------|
| `chr`, `pos`, `ref`, `alt`, `gene` | as above |
| `clnsig`     | Raw ClinVar CLNSIG (may be multi-valued). VEP encodes internal commas as `&`. |
| `clnrevstat` | Raw ClinVar CLNREVSTAT string (VEP `&`-encoded). |
| `condition`  | ClinVar disease name(s) from CLNDN, `; `-joined; `not_provided`/`not_specified` dropped when other terms exist. |
| `clnsigconf` | ClinVar CLNSIGCONF — the per-submitter breakdown when the call is conflicting (empty otherwise). |
| `stars`      | Star count derived from `clnrevstat` (0–4). |
| `is_clinvar_PLP` | 0/1. Requires P/LP CLNSIG, review-status ≥ `params.clinvar_min_stars`, and no benign co-classification. Conflicting calls are excluded (1 star). |

## `results/variants/acmg_plp.tsv`
Produced by **fastVEP** (Huang-lab/fastVEP) via its native `--acmg`
(Richards 2015 + ClinGen SVI), parsed by `acmg_fastvep/parse_fastvep_acmg.py`.
Per variant the most severe call across transcripts is kept (P > LP > VUS > LB > B).

| column | description |
|--------|-------------|
| `chr`, `pos`, `ref`, `alt`, `gene` | as above |
| `acmg_label` | fastVEP ACMG-AMP label: `Pathogenic` / `Likely_pathogenic` / `Uncertain_significance` / `Likely_benign` / `Benign`. |
| `acmg_criteria` | `;`-joined triggered ACMG codes, e.g. `PVS1;PM2_Supporting;PP3` (from fastVEP's `ACMG_CRITERIA` CSQ subfield). |
| `n_pathogenic_criteria` | Count of triggered pathogenic codes (PVS/PS/PM/PP). |
| `n_benign_criteria` | Count of triggered benign codes (BA/BS/BP). |
| `is_acmg_PLP` | 0/1 — 1 iff `acmg_label` is Pathogenic or Likely_pathogenic. |

## `results/variants/am_plp.tsv`
| column | description |
|--------|-------------|
| `chr`, `pos`, `ref`, `alt`, `gene` | as above |
| `am_score` | AlphaMissense pathogenicity score (0–1) or empty. |
| `min_strength` | Configured minimum PP3 evidence tier (`PP3_Supporting` / `PP3_Moderate` / `PP3_Strong` / `PP3_VeryStrong`). |
| `is_AM_PLP` | 0/1. |

## ACMG intermediates
The `acmg` classifier (fastVEP) writes a full annotated VCF
(`<chunk>.fastvep.vcf`, CSQ incl. `ACMG`/`ACMG_CRITERIA`) in the task work dir as
an intermediate; only the analysis-ready `variants/acmg_plp.tsv` is published. To
keep the full fastVEP VCF, run the standalone `acmg_fastvep/` tool (which retains
it) — see `acmg_fastvep/README.md`.

## `results/variants/qc_per_gene.tsv`
| column | description |
|--------|-------------|
| `gene` | VEP `SYMBOL`. |
| `n_annotated` | Number of annotated records in gene. |
| `n_clinvar_PLP`, `n_acmg_PLP`, `n_AM_PLP` | Per-framework P/LP counts. |

## `results/manifest.json`
See `assets/schemas/manifest.schema.json`. Keys:
`reference_build`, `vep.version`, `vep.cache`, `clinvar_release`,
`gnomad_version`, `alphamissense.data_version`,
`alphamissense.calibration_version`, `container_digests`, `pipeline_git_sha`.
