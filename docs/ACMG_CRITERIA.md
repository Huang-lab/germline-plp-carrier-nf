# ACMG-AMP criteria reference, and how this pipeline reports them

This is a reference for the `acmg_criteria` codes and other ACMG-related
columns produced by `acmg_fastvep/build_acmg_carriers.py` (and
`acmg_fastvep/parse_fastvep_acmg.py`). None of the classification logic lives
in this repo — every code and label here is computed by fastVEP
(Huang-lab/fastVEP) itself, using the standard published criteria. This
document only explains what those codes mean and how this pipeline reports
fastVEP's own output; it does not define or modify any classification rule.

## Where the classification actually happens

fastVEP's `--acmg` mode implements the ACMG-AMP framework from
[Richards et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/), with
[ClinGen SVI](https://clinicalgenome.org/working-groups/sequence-variant-interpretation/)
calibration (e.g. REVEL thresholds for PP3/BP4, PM2 downgraded to
"Supporting" strength). It evaluates all 28 criteria per transcript, using
real evidence (gnomAD population frequency, ClinVar, REVEL, SpliceAI, gene
constraint, and more), then combines the criteria that fire into one of five
tiers via the standard ACMG combining rules.

This pipeline's job is narrower: **retrieve** that already-computed result
for the correct transcript (MANE Select) and **match** it to sample
genotypes to find carriers. See `acmg_fastvep/build_acmg_carriers.py`'s
module docstring for exactly how the transcript is chosen.

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
