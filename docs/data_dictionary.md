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
| `is_acmg_PLP`   | 0/1    | 1 iff InterVar assigns Pathogenic or Likely pathogenic. |
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

## `results/annovar_intervar/` (ACMG intermediates)
Published only when the `acmg` classifier runs. Per chunk:
- `<chunk>.hg38_multianno.txt` — ANNOVAR `table_annovar.pl` multi-annotation
  output (refGene, clinvar, gnomAD, dbNSFP, avsnp).
- `<chunk>*.intervar` — raw InterVar classification (the ACMG evidence-code
  breakdown that `acmg_plp.tsv` is post-processed from).
These are provenance/debugging artifacts; the analysis-ready call is
`variants/acmg_plp.tsv`.

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
