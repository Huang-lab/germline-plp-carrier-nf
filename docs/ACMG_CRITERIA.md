# ACMG-AMP criteria reference, and how this pipeline reports them

This is a reference for the `acmg_criteria` codes and other ACMG-related
columns produced by `acmg_fastvep/build_acmg_carriers_pick.py` (the default
carrier builder), `acmg_fastvep/build_acmg_carriers.py` (the older builder that
reads MANE Select only) and `acmg_fastvep/parse_fastvep_acmg.py`. None of the
classification logic lives
in this repo — every code and label here is computed by fastVEP
(Huang-lab/fastVEP) itself, using the standard published criteria. This
document only explains what those codes mean and how this pipeline reports
fastVEP's own output; it does not define or modify any classification rule.

## Where the classification actually happens

fastVEP's `--acmg` mode implements the ACMG-AMP framework from
[Richards et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/), with
[ClinGen SVI](https://clinicalgenome.org/working-groups/sequence-variant-interpretation/)
calibration (e.g. REVEL thresholds for PP3/BP4, PM2 downgraded to
"Supporting" strength). It evaluates, per transcript, the criteria that the
supplementary-annotation databases it was given make possible (gnomAD
population frequency, ClinVar, REVEL, SpliceAI, gene constraint, and more),
then combines the criteria that fire into one of five tiers via the standard
ACMG combining rules.

Not every criterion is evaluated in a given run.
fastVEP's own `docs/ACMG.md` (v0.5.0) says "a default run reports 26", and the tables below list 28 codes.
In the JSON output of fastVEP v0.5.0 for one stop-gained variant, annotated with only gnomAD allele frequency, gnomAD gene constraint and an OMIM-format gene-disease table, 26 distinct codes were present (the 28 below except PP5 and BP6) and 11 of them had `evaluated: true` in at least one transcript: BP1, BP3, BP7, PM1, PM2_Supporting, PM3, PM4, PM5, PP2, PS1 and PVS1.
For the codes that were not evaluated, fastVEP's own explanation was that PS2 and PM6 need a trio VCF, PS3 needs curated functional evidence (`--functional-evidence`), PS4 needs case-control statistics, and PP1 needs a multi-generation pedigree.
`run_fastvep.sh --acmg` supplies none of these inputs.

This pipeline's job is narrower: **retrieve** that already-computed result
for the correct transcript and **match** it to sample genotypes to find
carriers.
The default builder, `build_acmg_carriers_pick.py`, takes the block that fastVEP
itself marked `PICK=1` for each (allele, gene).
The older `build_acmg_carriers.py` takes the MANE Select block only.
See each script's module docstring for exactly how the transcript is chosen.

## The five-tier label (`acmg_label`)

| Shorthand (fastVEP `ACMG` field) | `acmg_label` in our tables | Meaning |
|---|---|---|
| P | Pathogenic | |
| LP | Likely_pathogenic | |
| VUS | Uncertain_significance | |
| LB | Likely_benign | |
| B | Benign | |

`is_acmg_PLP` is 1 for P or LP, 0 otherwise.

## The criteria codes (`acmg_criteria`)

Each code may carry a strength suffix: `_Very_Strong`, `_Strong`,
`_Moderate`, `_Supporting` (no suffix = the code's default strength). Several
codes appear in `acmg_criteria` joined by `;`.

### Pathogenic evidence

| Code | Meaning |
|---|---|
| PVS1 | Null variant (nonsense, frameshift, canonical ±1/2 splice site, initiation codon, or single/multi-exon deletion) in a gene where loss-of-function is an established disease mechanism |
| PS1 | Same amino acid change as an established pathogenic variant |
| PS2 | De novo (confirmed) in a patient with the disease, no family history |
| PS3 | Well-established functional studies show a damaging effect |
| PS4 | Prevalence in affected individuals significantly increased over controls |
| PM1 | Located in a mutational hot spot / critical functional domain with no benign variation |
| PM2 | Absent, or at extremely low frequency, in population databases |
| PM3 | In trans with a pathogenic variant (recessive disorders) |
| PM4 | Protein length change (in-frame indel, stop-loss) in a non-repeat region |
| PM5 | Novel missense change at a residue where a different pathogenic missense change was previously seen |
| PM6 | Assumed de novo, without confirmed parentage |
| PP1 | Cosegregates with disease in family members |
| PP2 | Missense in a gene with a low rate of benign missense variation, where missense is a common disease mechanism |
| PP3 | Multiple computational/predictive evidence supports a deleterious effect (e.g. REVEL, SpliceAI) |
| PP4 | Patient's phenotype is highly specific for the gene |
| PP5 | Reputable source reports pathogenic, without independent evidence available (ClinGen SVI recommends against using this on its own) |

### Benign evidence

| Code | Meaning |
|---|---|
| BA1 | Allele frequency >5% in population databases (stand-alone benign) |
| BS1 | Allele frequency greater than expected for the disorder |
| BS2 | Observed in a healthy adult for a fully penetrant early-onset disorder |
| BS3 | Well-established functional studies show no damaging effect |
| BS4 | Lack of segregation in affected family members |
| BP1 | Missense in a gene where primarily truncating variants cause disease |
| BP2 | Observed in trans (recessive) or cis (dominant) with a pathogenic variant |
| BP3 | In-frame indel in a repetitive region with no known function |
| BP4 | Multiple computational/predictive evidence suggests no impact |
| BP5 | Found in a case with an alternate molecular explanation |
| BP6 | Reputable source reports benign, without independent evidence available (ClinGen SVI recommends against using this on its own) |
| BP7 | Synonymous variant with no predicted splice impact and no conservation signal |

## Other columns this pipeline adds

- **`consequence` / `impact` / `hgvsc` / `hgvsp` / `exon` / `intron`** — retrieved from
  the *same* MANE Select transcript block used for `acmg_label`/`acmg_criteria`,
  never mixed with a different transcript's context.
- **`n_pathogenic_criteria` / `n_benign_criteria`**: how many codes in `acmg_criteria` point towards pathogenic (PVS, PS, PM, PP) and towards benign (BA, BS, BP).
- **`is_lof_consequence`** — 1 if `consequence` is one of the standard fixed
  "predicted loss-of-function" terms (`stop_gained`, `frameshift_variant`,
  `splice_acceptor_variant`, `splice_donor_variant`, `start_lost` — the
  gnomAD/LOFTEE convention). This is a plain category-membership check, not
  a classification: it does **not** account for NMD escape or last-exon
  position, so it is the *basic* pLoF flag, not the more rigorous
  LOFTEE-style one. A variant can be `is_lof_consequence=1` and still be
  classified VUS/benign by fastVEP if the broader ACMG evidence does not
  support pathogenicity (gene mechanism, population frequency, etc.).
- **`spliceai_ds_ag/al/dg/dl`** — SpliceAI delta scores (acceptor gain/loss,
  donor gain/loss) for this allele, from fastVEP's separate `SpliceAI` INFO
  field (allele-keyed, not transcript-keyed).
- **`acmg_label = NO_MANE_TRANSCRIPT`** — this variant has no CSQ block with a
  MANE Select transcript at all (e.g. genuinely intergenic, or overlapping
  only non-coding/pseudogene transcripts, which MANE does not cover). The
  `gene` and `consequence` columns still list every gene symbol / consequence
  term seen across all of that variant's transcript blocks (semicolon-joined),
  so a flagged row can be told apart from truly intergenic space — nothing is
  picked or resolved on your behalf.
  This label is written by `build_acmg_carriers.py`.
  It is not written for a variant that lies near a MANE gene: such a variant is reported against that gene's MANE Select block, as `upstream_gene_variant` or `downstream_gene_variant` (VUS in the run described next).
  On a synthetic genome in which every gene was within 5 kb of another (fastVEP's default `--distance 5000`), four stop_gained variants in genes without a MANE Select transcript were therefore absent from the table, and the summary said `0 flagged NO_MANE_TRANSCRIPT`.
  `build_acmg_carriers_pick.py` keeps such variants.

## Per-gene pick variant (`build_acmg_carriers_pick.py`)

This is the default carrier builder.
It writes the same columns as above plus the ones listed below, from the same fastVEP output (no re-run).
Instead of MANE Select only, each (variant, gene) reports one CSQ block.
`--selection` says which:

- `fastvep_pick` (default): the block that fastVEP itself marked `PICK=1` for that (allele, gene).
  The input must be fastVEP v0.5.0 or later output made with `--flag-pick-allele-gene`, which `run_fastvep.sh --acmg` passes.
  Nothing is ported or re-implemented, and a CSQ without a `PICK` field is an error (exit 2).
- `python_pick` (DEPRECATED): a port of Huang-lab/fastVEP `crates/fastvep-annotate/src/pick.rs`, kept for a one-time comparison against output that has no `PICK` field (the v6 outputs).
  The tier table and `--pick-order` below apply to this mode only.
  `--pick-order` together with the default mode is an error.

### `python_pick` tiers

| tier | scores (lower wins) |
|---|---|
| `mane_select` / `mane_plus_clinical` | 0 if the field is set, else 1 |
| `canonical` | 0 if `CANONICAL=YES`, else 1 |
| `appris` | P1..P5 < principal < ALT1/2 < alternative < unrecognised < absent |
| `tsl` | the TSL number (1 best); absent = 255 |
| `biotype` | 0 if `protein_coding`, else 1 |
| `ccds` | 0 if a CCDS ID is set, else 1 |
| `rank` | most severe SO term's VEP 115 rank (transcript_ablation=1 … ); none = worst |

In `python_pick` mode the default is `--pick-order mane_select,mane_plus_clinical,canonical,tsl,biotype,ccds,rank`;
the VEP-path order `mane_select,canonical,tsl,biotype,rank,ccds` is also
accepted. Remaining ties go to the lowest transcript ID. Genes are grouped by
the CSQ `Gene` ID, as fastVEP does.

### Columns added by the pick table

- **`transcript_id`** — `Feature` of the picked block.
- **`pick_tier`**: `fastvep_pick` on every row in the default mode.
  In `python_pick` mode it is the first tier at which the winner beat the runner-up;
  `only_candidate` when the gene had one transcript block, `tie_id` when all
  listed tiers tied and the transcript ID decided.
- **`is_mane_select`** — 1 if the picked block is the MANE Select transcript.
- **`acmg_label_mane`**: the MANE Select block's label for that gene.
  If the gene has no MANE Select block at this variant it is `NO_MANE_TRANSCRIPT` in the default mode (the row still reports fastVEP's `PICK=1` block, so the call is kept), and empty in `python_pick` mode.
- **`mane_discordant`** — 1 if `is_acmg_PLP` differs between the picked block
  and the MANE block, 0 if it agrees, empty if there is no MANE block. With any
  order that starts with `mane_select`, a gene's MANE block always wins, so
  this is 0 or empty; it only becomes 1 with an order such as `rank,...`.
- With `--drop-flanking-genes` (default) a gene whose only consequences at the
  variant are `upstream_gene_variant`/`downstream_gene_variant` gets no row.
- **`acmg_label = NO_TRANSCRIPT`** — only for a variant with no gene-bearing
  transcript block at all (intergenic / regulatory only).
- Carrier `zygosity = hemi` — a haploid non-ref GT (e.g. `1` on male chrX
  outside the PARs). These were previously dropped; the fix applies to
  `build_acmg_carriers.py` too.
- **`acmg_disorders`**: the trailing `ACMG_DISORDERS` field of fastVEP v0.5.0 (`INHERITANCE:CLASS:CRITERIA` per mode of inheritance, joined by `&`), copied as it is.
- **`pick_flag`**: `ok` on a normal row.
  On a row that the PICK check flagged it says why: `no_pick`, or `multi_pick:` followed by each flagged transcript and its label.
- **`acmg_label = NO_PICK_FLAG` / `MULTIPLE_PICK_FLAGS`** (default mode): the (allele, gene) had no `PICK=1` block, or more than one.
  The row is not resolved, and `is_acmg_PLP` is 0.
  fastVEP's `--flag-pick-allele-gene` flags exactly one block per (allele, gene), so such a row means the output and this parser disagree.
  Both tables are still written, then the script exits 3 and writes a FAILED marker (a DONE marker otherwise).
  `--allow-pick-anomalies` accepts the rows, and `--no-markers` writes neither marker.
  In a test on a synthetic cohort, output made with `--flag-pick` or with `--flag-pick-allele` gave exit 3 with five `NO_PICK_FLAG` rows, and output made with `--pick` (which removes the `PICK` field) gave exit 2.

### Functional and population columns (pick variant only)

Appended after `mane_discordant`, so earlier column positions are unchanged.

From the **picked** transcript's CSQ block:
- **`protein_position`, `amino_acids`, `cds_position`** — e.g. `1443`, `R/*`, `4327`.
- **`flags`** — `cds_start_NF` / `cds_end_NF` mark a transcript whose coding
  sequence is incomplete; a pLoF call on one is a common false positive.
- **`existing_variation`** — known-variant IDs (rsIDs etc.), `&`-joined.
- **`is_last_exon`** — 1 if `exon` is the transcript's last exon (`N/N`), 0 if
  not, empty when there is no exon number. A stop in the last exon usually
  escapes nonsense-mediated decay. Only a partial NMD check: the "last 50 bp
  of the penultimate exon" part needs exon boundaries the CSQ does not carry.

From fastVEP's supplementary-annotation INFO fields, matched to this record's
alt allele (the VCF ALT, or VEP's anchor-trimmed spelling for indels, e.g.
`-` for a deletion). Records where such a field exists but has no entry for
the alt are counted in the stderr summary.
- **`clinvar_significance`, `clinvar_review_status`** — from `FV_CLINVAR`;
  several ClinVar records for one allele are `;`-joined.
- **`gnomad_af`** — `ALL_AF`; **`gnomad_hom_count`** — `ALL_HC`.
- **`gnomad_grpmax_af`, `gnomad_grpmax_pop`** — highest AF among AFR, AMR,
  EAS, NFE and SAS: gnomAD's "grpmax", which leaves out the bottlenecked
  groups (Ashkenazi Jewish, Finnish, Middle Eastern, Amish, remaining). This
  is the frequency BA1/BS1 refer to.
- **`revel`** — `FV_REVEL` score (the maximum if several).
- **`gene_pli`, `gene_loeuf`** — gnomAD constraint for this row's gene (matched
  by symbol) from `FV_GNOMAD_GENE`. Low LOEUF (< ~0.6) = loss-of-function
  intolerant.
