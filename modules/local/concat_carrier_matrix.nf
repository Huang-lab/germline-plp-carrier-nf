process CONCAT_CARRIER_MATRIX {
    input:
    path(chunk_tsvs)

    output:
    path("carrier_matrix.tsv"), emit: long_tsv

    script:
    """
    set -euo pipefail
    # Write header from the first file; concat body from all files.
    files=( ${chunk_tsvs} )
    head -n1 "\${files[0]}" > carrier_matrix.tsv
    for f in "\${files[@]}"; do
        tail -n +2 "\$f" >> carrier_matrix.tsv
    done
    """
}
