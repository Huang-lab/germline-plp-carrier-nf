process ACMG_FASTVEP {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf)
    path gff3
    path reference_fasta
    path sa_dir

    output:
    path("${chunk_id}.acmg_plp.tsv"), emit: tsv

    script:
    def allow_stub = (params.allow_stub_acmg as boolean)
    def fastvep    = params.fastvep_bin ?: 'fastvep'
    """
    set -euo pipefail
    if command -v ${fastvep} >/dev/null 2>&1 && [ -s "${gff3}" ] && [ -d "${sa_dir}" ]; then
        # ACMG-AMP via fastVEP (Huang-lab/fastVEP): its native --acmg runs
        # Richards 2015 + ClinGen SVI. fastVEP re-annotates from the variant
        # records (it ignores/strips any existing CSQ), so feeding the pipeline's
        # already-annotated VCF is fine. Needs a GFF3 + FASTA + supplementary
        # DBs (--sa-dir): gnomAD/ClinVar/REVEL/gene-level (fastvep sa-build).
        FASTA_ARG=""
        [ -s "${reference_fasta}" ] && FASTA_ARG="--fasta ${reference_fasta}"
        ${fastvep} annotate \\
            --input ${vcf} --output ${chunk_id}.fastvep.vcf \\
            --gff3 ${gff3} \$FASTA_ARG --sa-dir ${sa_dir} \\
            --acmg --output-format vcf
        # Collapse per-transcript ACMG to one row/variant (most severe) → pipeline schema.
        python3 ${projectDir}/acmg_fastvep/parse_fastvep_acmg.py \\
            --vcf ${chunk_id}.fastvep.vcf --out ${chunk_id}.acmg_plp.tsv
    elif [ "${allow_stub}" = "true" ]; then
        # Synthetic stub — ONLY for the test profile (params.allow_stub_acmg=true).
        # Emits a valid acmg_plp.tsv so DAG wiring/publishing is exercised without
        # a fastVEP install or supplementary databases.
        python3 ${projectDir}/bin/synthetic_acmg_fastvep.py --in ${vcf} --out ${chunk_id}.acmg_plp.tsv
    else
        echo "ERROR: ACMG requested but fastVEP is not available/configured." >&2
        echo "  Install fastvep (acmg_fastvep/setup_fastvep.sh) and build supplementary" >&2
        echo "  databases, then set params fastvep_gff3, fastvep_sa_dir (and optionally" >&2
        echo "  fastvep_bin / container_fastvep). See acmg_fastvep/README.md." >&2
        echo "  (For pipeline testing only, set --allow_stub_acmg true to emit a synthetic stub.)" >&2
        exit 1
    fi
    """
}
