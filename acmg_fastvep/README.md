# acmg_fastvep — standalone fastVEP v0.5.0 ACMG on QC'd VCFs (runbook)

Run **[fastVEP](https://github.com/Huang-lab/fastVEP)** v0.5.0 with ACMG-AMP
classification on the VCFs the Nextflow pipeline QC-normalizes
(`results-chrN/norm_qc/<chunk>.norm.vcf.gz`), **outside** `main.nf`, as plain
scripts and LSF jobs. fastVEP is a native binary: no container, no `module load`
for fastVEP itself (samtools/htslib are needed for setup and fetching only).

Everything is pinned. One versioned directory holds the binary, the fastVEP
source at the tag (for its helper scripts), the gene model, the transcript cache,
the downloaded sources and the built databases; every job reads it through one
`env.sh`, and every output records the hashes it was made with.

```
$FV_HOME   (default /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/variants_PLP/data/fastvep-v0.5.0)
├── bin/fastvep                 v0.5.0 release binary (sha256 in env.sh / fastvep.provenance)
├── src/                        fastVEP @ v0.5.0 (19a159f9): scripts/, analysis/.../sa_sources/
├── ref/                        Homo_sapiens.GRCh38.115.gff3 + .fastvep.cache
├── sa_src/                     downloads + SOURCES.tsv (url, date, size, sha256)
├── sa_db/                      the 100 databases, each with .done; SA_DB_SHA256.txt, SA_DB_MANIFEST.tsv
├── fastvep.provenance          tag, commit, tarball + binary sha256, --version, glibc, host
└── env.sh                      FV_HOME FV_BIN FV_BIN_SHA256 FV_GFF3 FV_FASTA FV_CACHE FV_SA FV_SRC
```

The old install (fastvep 0.3.0 on PATH at
`/sc/arion/projects/rg_huangk06/variants_PLP_MSM/tools/bin/fastvep`, databases under
`.../variants_PLP/data/fastvep`) is not touched and not used. The scripts never run
a `fastvep` found on PATH: they run `$FV_BIN` and check its sha256.

## Files

| File | Runs on | What |
|---|---|---|
| `setup_fastvep.sh` | login | install binary + source tree, GFF3, FASTA index; cache via bsub; write `env.sh` |
| `fetch_sa_sources.sh` | login | download and verify every source into `sa_src/` |
| `build_sa_dbs.sh` / `build_sa_dbs.lsf` | compute | build `sa_db/` (`--single` job + `--chrom` array) |
| `check_sa_db.sh` | compute (full) / anywhere (`--quick`) | exact file set, `check_acmg_stack.py`, seal with SA_DB_SHA256.txt |
| `run_fastvep.sh` | compute | annotate one chunk; `--acmg` adds `--flag-pick-allele-gene` |
| `run_fastvep_batch.sh` / `run_fastvep.lsf` | compute | all chunks of a `norm_qc/` dir, with provenance |
| `parse_fastvep_acmg.py` | (called) | CSQ → `<chunk>.acmg_plp.tsv` (MOST SEVERE call over ALL transcripts; not a carrier source, see below) |
| `fv_lock.sh` | (sourced) | one writer per chunk / per batch outdir (flock, else mkdir lock) |
| `build_acmg_carriers_pick.py` | (called) | per (allele, gene) row + carriers; default `--selection fastvep_pick` |
| `run_acmg_carriers_chr.sh` | compute | carriers for one chromosome dir (default `--mode fastvep_pick`) |
| `compare_acmg_pick_vs_mane.py` | anywhere | MANE-only vs pick results |
| `extract_gene_carriers.py` | anywhere | carriers of a gene list across chromosomes |
| `build_acmg_carriers.py` | (called) | MANE Select only (`--mode mane`), for v6 comparisons |

## 0. Setup (login node, once)

```bash
cd germline-plp-carrier-nf            # this repo / worktree
acmg_fastvep/setup_fastvep.sh \
    --fasta /sc/arion/projects/rg_huangk06/variants_PLP_MSM/refs/vep_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa
#   [--fv-home <dir>]  [--from-source]  [--allow-path-mismatch]
```

* Downloads `fastVEP-v0.5.0-x86_64-unknown-linux-gnu.tar.gz`, checks its sha256
  (`3b081b9e…ce82`), installs `bin/fastvep` (built for glibc ≥ 2.34; Minerva has
  2.34). `--from-source` builds the tag with `cargo build --release --locked`
  instead. Either way `src/` is cloned at the tag and its commit asserted.
* If a different `fastvep` is first on PATH (the 0.3.0 one is), setup **stops**
  and says so. Take it off PATH, or pass `--allow-path-mismatch` once you have
  read why (jobs use `$FV_BIN`; only hand-typed `fastvep` would get the old one).
* The transcript cache needs ~16 GB, so setup prints a `bsub ... --cache-only`
  command instead of building it on the login node. Submit it, wait for it.
  (`--build-cache` builds it in place, e.g. in an interactive compute session.)

Then, in every shell and job: `source $FV_HOME/env.sh` (jobs: `export FV_ENV=$FV_HOME/env.sh`).

## 1. Fetch the sources (login node)

```bash
source $FV_HOME/env.sh
acmg_fastvep/fetch_sa_sources.sh --clinvar-date 20261109    # pin a weekly ClinVar; omit for the latest
```

ClinVar VCF (md5 + `gzip -t` + AF_EXAC header), `variant_summary.txt.gz`,
**gnomAD v4.1 exomes** chr1–22, X, Y (~198 GB; BGZF EOF block + `tabix -l`;
`--full-check` adds `bgzip -t`), constraint metrics, REVEL v1.3 (split per
chromosome), **ClinGen Gene-Disease Validity** (dated by download day), UCSC rmsk.
Resumable; re-run after an interruption. `--skip-gnomad`, `--chroms "21 22"`.
OMIM is opt-in only: `--omim-genemap2 /path/licensed/genemap2.txt`.
Result: `sa_src/SOURCES.tsv`. `variant_summary.txt.gz` (for `clinvar_protein`) is
NCBI's **moving** file with no weekly copy: the copy in `sa_src` is the pinned one,
its server date goes to `sa_src/VARIANT_SUMMARY_DATE` next to `CLINVAR_DATE`, and
its sha256 is in SOURCES.tsv and in `clinvar_protein.oga.done`. Refreshing ClinVar
means refreshing both. A re-run never re-downloads a file that already passes its
check, and a resume that the server answers with HTTP 416 counts as complete.

## 2. Build the databases (compute: one array + one job)

```bash
export MINERVA_ALLOCATION=acc_<project> FV_ENV=$FV_HOME/env.sh
mkdir -p logs
# per chromosome: index 1-22 -> chr1-22, 23 -> X, 24 -> Y. 1 core, 8 GB, 4 h each.
bsub -P "$MINERVA_ALLOCATION" < acmg_fastvep/build_sa_dbs.lsf
# whole-genome databases: 1 core, 16 GB, 2 h.
FV_SA_MODE=single bsub -P "$MINERVA_ALLOCATION" -J fv050-sa-single \
    -R "rusage[mem=16000] span[hosts=1]" -W 2:00 \
    -oo logs/fv050-sa-single.%J.out -eo logs/fv050-sa-single.%J.err \
    < acmg_fastvep/build_sa_dbs.lsf
```

Exactly as v0.5.0 `docs/ACMG_SETUP.md`: `gnomad_chrN`, `spliceai_chrN` and
`phylop_chrN` (SpliceAI and PhyloP distilled from the same gnomAD file by
`scripts/extract_gnomad_scores.py`) with `--chunk-bits 16`; `revel_chrN` (no Y:
REVEL has none); `clinvar`; `clinvar_protein` from `variant_summary`;
`gnomad_genes`; `omim` from ClinGen via `clingen_gdv_to_oga.py` (`FV_SA_ARGS=--use-omim`
for licensed OMIM); `repeatmasker` (`--source custom_bed --name repeatmasker`).
`sa-build --format auto` writes `.osa2` for every allele-level source (clinvar and
revel included), `.oga` gene-level, `.osi` for RepeatMasker.

A database reaches `sa_db/` only after `sa-build` exits 0, then gets `<db>.done`,
which records the upstream source file and its sha256, the binary sha256 and the
sa-build args. Re-submitting skips a database only when all of those still match;
a new ClinGen export, a re-split REVEL, a new binary or `--use-omim` rebuilds the
affected databases without `--force` (stderr says why). Never `--allow-empty`. Logs:
`$FV_HOME/sa_build_logs/`. To reclaim the 198 GB once a chromosome is built:
`FV_SA_ARGS=--delete-gnomad-src` on the array (deletes only after gnomad, spliceai
and phylop of that chromosome all have `.done`).

**Defaults that matter:** ClinGen GDV, not OMIM, behind the `omim` key. gnomAD
v4.1 **exomes**. **Never put joint or genomes gnomAD databases in the same
sa_db**: fastVEP keeps one `gnomad` answer per allele (the database sorting
last), so mixing silently swaps frequencies. `check_sa_db.sh` refuses it.

## 3. Check and seal (compute, once per build)

```bash
export FV_ENV=$FV_HOME/env.sh        # check_sa_db.sh sources it; bsub passes the env on
bsub -P "$MINERVA_ALLOCATION" -J fv050-sa-check -n 1 -R "rusage[mem=16000]" -W 1:00 \
    -oo logs/fv050-sa-check.%J.out -eo logs/fv050-sa-check.%J.err \
    acmg_fastvep/check_sa_db.sh
```

Fails unless `sa_db/` holds exactly the 100 expected databases, each with
`.done`, and nothing else a fastVEP run would load; every gnomAD database came from
v4.1 exomes; v0.5.0 `scripts/check_acmg_stack.py` sees all nine sources answer;
and the ClinGen source knows BRCA1, BRCA2, TP53, PTEN, MYH7. Only then it writes
`SA_DB_SHA256.txt` + `SA_DB_MANIFEST.tsv`. Every run checks `--quick` against
that manifest first, and records `sha256(SA_DB_SHA256.txt)` as its `sa_db_sha256`.

## 4. chr21 smoke test

```bash
export FV_ENV=$FV_HOME/env.sh FV_MODE=acmg
export FV_IN=<batch>/results-chr21/norm_qc  FV_OUT=<batch>/results-fastvep-chr21_fv050
bsub -P "$MINERVA_ALLOCATION" -J fv050-acmg-chr21 < acmg_fastvep/run_fastvep.lsf
```

What to verify before anything else runs:

1. **No gated warnings.** fastVEP exits 0 on all of these, so a chunk whose stderr
   says `incomplete reference stack`, `was not loaded`, `hold nothing of
   chromosomes`, `answered by more than one database`, `could not be read for`
   (an SA lookup I/O error: no annotation, PM2 can fire), `could not be read and
   were not annotated` (records dropped), `not a contig of the gene models` or
   `are on contigs that the gene models and the reference do not have` (records
   called intergenic) fails and publishes nothing. So: every line of
   `$FV_OUT/batch_manifest.tsv` is `ok`, and
   `grep -lE 'reference stack|was not loaded|hold nothing|more than one database|could not be read|not a contig|are on contigs' $FV_OUT/*.fastvep.stderr.log`
   prints nothing.
2. **CSQ fields.** `zcat $FV_OUT/<chunk>.fastvep.vcf.gz | grep -m1 '^##INFO=<ID=CSQ' | sed 's/.*Format: //; s/">$//' | tr '|' '\n' | grep -xcE 'MANE_SELECT|ACMG|ACMG_CRITERIA|ACMG_DISORDERS|PICK'`
   prints `5` (ACMG_DISORDERS is the last field, PICK follows FLAGS).
3. **Exactly one PICK=1 per allele + gene.** Run the carriers step (§5); its
   `summary.txt` must say `with no PICK=1: 0; with more than one: 0`. Anything
   else is listed on stderr (`chunks.stderr.log`), those rows are flagged
   `NO_PICK_FLAG` / `MULTIPLE_PICK_FLAGS` (never resolved by a pick of our own),
   and the carriers job **fails (exit 5) after writing its tables**: per fastVEP's
   pick.rs it cannot happen, so a non-zero count means parser and output disagree.
   `--allow-pick-anomalies` accepts it once understood.
4. **Provenance.** `$FV_OUT/fastvep_run_provenance.txt` names the binary sha256 and
   version (0.5.0), cache sha256, FASTA, `sa_db_sha256`, the pipeline git rev and
   LSB_JOBID; each chunk has `<chunk>.fastvep.provenance` with the same hashes.

## 5. Carriers and the NO_MANE / pick comparison on chr21

```bash
acmg_fastvep/run_acmg_carriers_chr.sh --mode fastvep_pick -j 4 \
    --in-dir <batch>/results-fastvep-chr21_fv050 \
    -o <batch>/results-carriers-fastvep_pick-chr21_fv050
cat <batch>/results-carriers-fastvep_pick-chr21_fv050/summary.txt
```

Each (variant, gene) row is the block fastVEP itself flagged `PICK=1` (its
default pick order: mane_select, mane_plus_clinical, canonical, appris, tsl,
biotype, ccds, rank); `acmg_label_mane` is the MANE Select block's call or
`NO_MANE_TRANSCRIPT`; `mane_discordant=1` where the two disagree on P/LP;
`acmg_disorders` carries v0.5.0's per-inheritance-mode calls.

Compare on chr21, before the full run:

```bash
# fv050 fastvep_pick vs the v6 MANE-only result (different prefixes and suffixes)
python3 acmg_fastvep/compare_acmg_pick_vs_mane.py \
    --mane-root <v6 MANE results root> --pick-root <batch> \
    --pick-prefix results-carriers-fastvep_pick- --strip-dir-suffix \
    --out <batch>/comparison-fv050-chr21
# how many rows have no MANE transcript, and how many of them are P/LP
awk -F'\t' 'NR>1 && $25=="NO_MANE_TRANSCRIPT" {n++; p+=($10==1)} END {print n+0" rows, "p+0" P/LP"}' \
    <batch>/results-carriers-fastvep_pick-chr21_fv050/acmg_variants.tsv
# one-time: does the DEPRECATED Python pick port choose the same transcripts on the same output?
acmg_fastvep/run_acmg_carriers_chr.sh --mode pick -j 4 \
    --in-dir <batch>/results-fastvep-chr21_fv050 -o <scratch>/results-carriers-pick-chr21_fv050
```

`pick_tier_counts.tsv`, `mane_discordant_plp.tsv` and `plp_only_in_one.tsv` in
the comparison dir are what to read. **Do not change the thresholds in
`docs/ACMG_CRITERIA.md` until the chr21 comparison has been reviewed.**

## 6. Full run

From the MSM repo: `minerva_jobs/submit_fastvep_fv050_all_chr.sh` submits, per
chromosome, `run_fastvep.lsf` → `results-fastvep-chrN_fv050/` and, when it is
done, the carriers job → `results-carriers-fastvep_pick-chrN_fv050/` (job names
`fv050-*`). Or all chromosomes' carriers as one array:

```bash
acmg_fastvep/submit_acmg_carriers.sh --batch-root <batch> --dir-suffix _fv050 \
    --mode fastvep_pick --out-root <batch>/results-acmg-fv050
```

Re-running the batch skips a finished chunk only if every key field of its
provenance matches: mode, binary sha256, sa_db hash, input chunk sha256, transcript
cache sha256, FASTA path, acmg-config sha256 (or none) and the fastVEP flags.
Output made with anything else is refused, naming the fields that differ
(`--force` / `FV_FORCE=1` redoes it). Use a fresh output dir per sa_db build.

Only one job writes a chunk at a time: partial files are named per job
(`.<chunk>.*.partial.<LSB_JOBID>.<pid>`), each chunk takes `<outdir>/.<chunk>.lock`
and a local batch takes `<outdir>/.batch.lock` (flock on Minerva; a mkdir lock
elsewhere, whose error names the holder and how to remove a stale one). A second
writer exits 7.

### Gene panel

The panel comes from the MSM repo, `scripts/gene_panel.py` (`GENE_PANEL`, the
single source of truth). Do not hard-code it here or in job scripts:

```bash
python3 -c 'import sys; sys.path.insert(0, "<MSM repo>/scripts"); from gene_panel import GENE_PANEL; print("\n".join(GENE_PANEL))' > panel_genes.txt
python3 acmg_fastvep/extract_gene_carriers.py --pick-root <batch> \
    --prefix results-carriers-fastvep_pick- --genes panel_genes.txt \
    --out panel_carriers.tsv --out-summary panel_summary.tsv
```

## Plain annotation (no ACMG)

`run_fastvep.sh -i <chunk> -o <dir>` without `--acmg` needs only env.sh (binary,
cache, FASTA): `<chunk>.fastvep.vcf.gz`. `--pick` is allowed here only (it and
`--flag-pick-allele-gene` are one exclusive option group in fastVEP).

## Notes

* **Contigs:** the pipeline's `*.norm.vcf.gz` use Ensembl names (`1,2,…`), as do the
  GFF3, cache and FASTA; fastVEP resolves the `chr`-prefixed gnomAD/UCSC sources.
* `run_fastvep.sh` passes `--transcript-cache` and not `--gff3` (fastVEP ignores the
  GFF3 when an explicit cache loads); FV_GFF3 is recorded in provenance only.
* **`<chunk>.acmg_plp.tsv` and `all.acmg_plp.tsv` (`--concat`) are NOT a carrier
  source.** `parse_fastvep_acmg.py` keeps, per variant, the most severe ACMG call
  over ALL transcripts and genes (the max-over-transcripts behaviour the pick work
  replaced), so a variant is P/LP there if any transcript says so. Carriers come
  from `run_acmg_carriers_chr.sh --mode fastvep_pick` (`acmg_variants.tsv` /
  `acmg_carriers.tsv`). The TSV keeps the pipeline schema (`chr, pos, ref, alt,
  gene, acmg_label, acmg_criteria, n_pathogenic_criteria, n_benign_criteria,
  is_acmg_PLP`) for comparison only; do not feed it to `build_carrier_matrix.sh`
  for carrier counts.
* v6 outputs (fastVEP pre-0.5.0, no PICK field) still work with `--mode mane` and
  the deprecated `--mode pick`; `--mode fastvep_pick` refuses them.

## Tests

```bash
python -m pytest acmg_fastvep/tests/ -q
```

fastvep and the fastVEP helper scripts are stubbed in every test except
`test_run_fastvep_smoke.py`, which runs only with `FV_ENV` set (on Minerva). The
carriers driver tests need GNU `realpath` (`-e`/`-m`).
