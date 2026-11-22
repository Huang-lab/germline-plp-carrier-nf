# ACMG carrier table of record

Issue #9 found that nothing in this repository said which carrier table the
lab's downstream analyses (penetrance, panel counts) read, or which code and
mode made it. This file is that record. Change it, in a commit, whenever the
table of record changes.

## Decision (2026-11-22)

The ACMG carrier table of record is the output of

```bash
acmg_fastvep/run_acmg_carriers_chr.sh --mode fastvep_pick ...
```

(`--mode fastvep_pick` is the default), run on **fastVEP v0.5.0** output
(tag `v0.5.0`, commit `19a159f94fdb3d2aeac51812209387509e26814a`, pinned by
`setup_fastvep.sh`), with the pipeline code on `main` at
`d1d83ee` (the merge of `fastvep-v0.5.0`, PR #14) or later.

Per chromosome, the table is the pair `acmg_variants.tsv` / `acmg_carriers.tsv`
in `results-carriers-fastvep_pick-chrN_fv050/` (README sections 5 and 6).

What `fastvep_pick` means: for each (variant, gene), the CSQ block that
fastVEP itself flagged `PICK=1`, in fastVEP's pick order (mane_select,
mane_plus_clinical, canonical, appris, tsl, biotype, ccds, rank). For a gene
with a MANE Select transcript that is the MANE Select call.
`docs/ACMG_CRITERIA.md` lists what that call evaluates.

## Not the table of record

| Output | Why not |
|---|---|
| `<chunk>.acmg_plp.tsv`, `all.acmg_plp.tsv` | most severe call over **all** transcripts and genes; for comparison only (README, Notes) |
| `--mode mane` | v6 (pre-0.5.0) outputs only; MANE-only, so genes without MANE Select get `NO_MANE_TRANSCRIPT` |
| `--mode pick` | deprecated Python port of the pick; kept for a one-time comparison |
| any table made with fastVEP before v0.5.0 (the v6 run) | different binary and annotation databases |

## How a given table proves it is one

Every carrier directory written by `run_acmg_carriers_chr.sh` holds either
`FAILED` (do not use) or `DONE`. Since this change, `DONE`'s `detail` line
records:

* `mode=fastvep_pick`
* `pipeline_rev=<40-char commit>` and `dirty=no` (with `dirty=yes` the commit
  alone does not reproduce the table)
* `chunks=<n> inventory=checked` (every chunk the annotation batch wrote is
  present, none failed, all from one binary and database; PR #8)

The annotation side is recorded next to the input chunks:
`fastvep_run_provenance.txt` (binary sha256, sa_db hash, pipeline rev, args)
and one `<chunk>.fastvep.provenance` per chunk.

A table whose `DONE` lacks `pipeline_rev` was made before this change, and its
commit is not recorded.

## The cohort run behind the numbers quoted

Fill this in once the full v0.5.0 run is done, before any ACMG number is quoted
(issue #9, suggested fix 3):

| Field | Value |
|---|---|
| batch root | _(path on Minerva)_ |
| carrier dirs | `results-carriers-fastvep_pick-chr{1..22,X,Y}_fv050/` |
| `pipeline_rev` / `dirty` from `DONE` | _(fill in)_ |
| `sa_db_sha256` from `fastvep_run_provenance.txt` | _(fill in)_ |
| date of the run | _(fill in)_ |
| checked against the v6 MANE table with `compare_acmg_pick_vs_mane.py` | _(yes/no, where the comparison is)_ |

## Annotation-database notes for this decision

* **SpliceAI covers indels.** `build_sa_dbs.sh` builds `spliceai_chrN` from
  gnomAD v4.1 exomes' own `spliceai_ds_max` with fastVEP v0.5.0
  `scripts/extract_gnomad_scores.py`, which writes every ALT allele that has a
  score, with no SNV-only filter. gnomAD v4.1 exomes carries the score on
  indels: in three samples of the chr21 file, 2,964 of 5,208 indel records had
  one (27,318 of 59,361 SNVs did). So only indels gnomAD has seen get a score;
  the separate indel merge in PR #5 (`build_spliceai_indels.sh`, for the v6
  database) is not needed for this table.
* `sa_db/SA_DB_SHA256.txt` is written by `check_sa_db.sh` for the v0.5.0 build
  (100 databases) and required by `run_fastvep.sh`. It is not the v6
  checksum list in PR #5.
