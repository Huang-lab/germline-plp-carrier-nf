process CLINVAR_CLASSIFY {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf)

    output:
    path("${chunk_id}.clinvar_plp.tsv"), emit: tsv

    script:
    """
    set -euo pipefail
    ${projectDir}/bin/classify_clinvar.py \\
        --vcf ${vcf} \\
        --min-stars ${params.clinvar_min_stars} \\
        --out ${chunk_id}.clinvar_plp.tsv
    """
}
