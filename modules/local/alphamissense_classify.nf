process ALPHAMISSENSE_CLASSIFY {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf)
    path am_calibration_tsv

    output:
    path("${chunk_id}.am_plp.tsv"), emit: tsv

    script:
    """
    set -euo pipefail
    ${projectDir}/bin/classify_alphamissense.py \\
        --vcf ${vcf} \\
        --calibration ${am_calibration_tsv} \\
        --min-strength ${params.am_min_strength} \\
        --default-threshold ${params.am_default_threshold} \\
        --out ${chunk_id}.am_plp.tsv
    """
}
