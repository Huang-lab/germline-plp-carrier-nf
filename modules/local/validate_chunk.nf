process VALIDATE_CHUNK {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf)
    val classifiers_csv

    output:
    tuple val(chunk_id), path(vcf), emit: vcf
    tuple val(chunk_id), path("${chunk_id}.validate.json"), emit: report

    script:
    """
    set -euo pipefail
    ${projectDir}/bin/validate_chunk.py \\
        --vcf ${vcf} \\
        --min-symbol-fraction ${params.min_symbol_fraction} \\
        --classifiers "${classifiers_csv}" \\
        > ${chunk_id}.validate.json
    """
}
