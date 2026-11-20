# Real fastVEP v0.5.0 output for a synthetic cohort

`real_fastvep_v0.5.0_cohort.vcf` and `real_fastvep_v0.5.0_multiallelic_split.vcf` are what the released fastVEP v0.5.0 binary wrote for an invented genome and ten invented samples.
There is no participant data and no real gene model in them.
They are used by `acmg_fastvep/tests/test_real_fastvep_output.py`, which compares the carrier builders with a hand count written in that file.
The VCF files are force-added: the repository ignores `*.vcf`, and `hooks/pre-commit` exempts only `tests/synthetic/`, which is why they live here.

`make_synthetic_cohort.py` writes the genome (`synth.fa`), the gene models (`synth.gff3`), the input cohort (`cohort.vcf`), a truth table and three tiny supplementary-annotation sources.
It is seeded and writes the same files on every run.

To regenerate the fixtures (macOS arm64 shown; the release asset name changes with the platform):

```bash
python3 make_synthetic_cohort.py work && cd work && samtools faidx synth.fa
curl -sSLO https://github.com/Huang-lab/fastVEP/releases/download/v0.5.0/fastVEP-v0.5.0-aarch64-apple-darwin.tar.gz
shasum -a 256 fastVEP-v0.5.0-aarch64-apple-darwin.tar.gz   # bc1480656f33448ac784dc169b3928caa38ed70e3381384ac57d6596ca290bb3
tar xzf fastVEP-v0.5.0-aarch64-apple-darwin.tar.gz && FV=$PWD/fastVEP-v0.5.0-aarch64-apple-darwin/fastvep
mkdir sa_db
$FV sa-build --source omim        -i sa_src/genemap2.txt   -o sa_db/omim        --assembly GRCh38
$FV sa-build --source gnomad_genes -i sa_src/constraint.tsv -o sa_db/gnomad_genes --assembly GRCh38
$FV sa-build --source gnomad      -i sa_src/gnomad.vcf     -o sa_db/gnomad      --assembly GRCh38
$FV cache --gff3 synth.gff3 --fasta synth.fa -o synth.cache
# the flags acmg_fastvep/run_fastvep.sh --acmg uses:
$FV annotate --input cohort.vcf --output real_fastvep_v0.5.0_cohort.vcf --output-format vcf \
    --transcript-cache synth.cache --fasta synth.fa --no-progress --sa-dir sa_db --acmg --flag-pick-allele-gene
# the split form of the multi-allelic record, made by the pipeline's NORM_QC call:
grep -E '^#|^1\t16054\t' cohort.vcf > multi.vcf
bcftools norm -m -any -f synth.fa multi.vcf -Ov -o multi_split.vcf
$FV annotate --input multi_split.vcf --output real_fastvep_v0.5.0_multiallelic_split.vcf --output-format vcf \
    --transcript-cache synth.cache --fasta synth.fa --no-progress --sa-dir sa_db --acmg --flag-pick-allele-gene
```

The genome is deliberately gene-dense (genes less than 5 kb apart), so most variants also get an `upstream_gene_variant` or `downstream_gene_variant` block for a neighbouring gene; that is what the flanking-drop path is tested on.
