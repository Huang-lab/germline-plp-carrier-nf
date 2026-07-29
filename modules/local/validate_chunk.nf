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
    # Validation is a QC REPORT, not a gate: a low SYMBOL fraction is expected in
    # gene-sparse regions (e.g. chrX/chrY deserts, PAR boundaries) and must not
    # kill the run. The script still writes ok:true/false into the JSON, so bad
    # chunks are recorded and reviewable — we just don't let a false stop the
    # pipeline. `|| true` swallows the non-zero exit; the JSON is always written.
    ${projectDir}/bin/validate_chunk.py \\
        --vcf ${vcf} \\
        --min-symbol-fraction ${params.min_symbol_fraction} \\
        --classifiers "${classifiers_csv}" \\
        > ${chunk_id}.validate.json || true

    # Guarantee the report exists even if the script died before writing (so the
    # downstream join never breaks on a missing file).
    if [ ! -s ${chunk_id}.validate.json ]; then
        echo '{"ok": false, "errors": ["validate_chunk.py produced no output"], "n_records": 0}' \\
            > ${chunk_id}.validate.json
    fi
    """
}
