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
| `is_clinvar_PLP`| 0/1    | 1 iff variant is P/LP per ClinVar with review-status ≥ `params.clinvar_min_stars`. |
| `is_acmg_PLP`   | 0/1    | 1 iff InterVar assigns Pathogenic or Likely pathogenic. |
| `is_AM_PLP`     | 0/1    | 1 iff AlphaMissense score meets the gene-specific threshold at `params.am_min_strength`. |

Rows exist only for carriers (at least one alt allele in the person's genotype)
whose variant is qualifying under at least one of the three definitions.

## `results/carriers/carrier_matrix.wide.tsv` (optional pivot)
Rows: `variant_key = chr:pos:ref:alt`. Columns: samples. Cell = 1 iff carrier
of a variant that is P/LP under any framework, else 0.

## `results/variants/clinvar_plp.tsv`
| column | description |
|--------|-------------|
| `chr`, `pos`, `ref`, `alt`, `gene` | as above |
| `clnsig`     | Raw ClinVar CLNSIG (may be multi-valued, comma/underscore joined). |
| `clnrevstat` | Raw ClinVar CLNREVSTAT string. |
| `stars`      | Star count derived from `clnrevstat` (0–4). |
| `is_clinvar_PLP` | 0/1. |

## `results/variants/acmg_plp.tsv`
| column | description |
|--------|-------------|
| `chr`, `pos`, `ref`, `alt`, `gene` | as above |
| `acmg_label` | InterVar label: `Pathogenic` / `Likely pathogenic` / `Uncertain significance` / `Likely benign` / `Benign`. |
| `is_acmg_PLP` | 0/1. |

## `results/variants/am_plp.tsv`
| column | description |
|--------|-------------|
| `chr`, `pos`, `ref`, `alt`, `gene` | as above |
| `am_score` | AlphaMissense pathogenicity score (0–1) or empty. |
| `min_strength` | Configured minimum PP3 evidence tier (`PP3_Supporting` / `PP3_Moderate` / `PP3_Strong` / `PP3_VeryStrong`). |
| `is_AM_PLP` | 0/1. |

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
