process CARRIER_MATRIX {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(clinvar_tsv), path(acmg_tsv), path(am_tsv), path(gt)
    path keep_samples
    val classifiers_csv

    output:
    path("${chunk_id}.carrier_matrix.tsv"), emit: long_tsv

    script:
    // GT already extracted by CARRIER_GT (bcftools C-speed). Here we only join
    // genotypes to classifications and derive zygosity — pure Python, no bcftools.
    def has_cv = clinvar_tsv?.size() > 0 ? "--clinvar ${clinvar_tsv}" : ""
    def has_ac = acmg_tsv?.size()    > 0 ? "--acmg ${acmg_tsv}"       : ""
    def has_am = am_tsv?.size()      > 0 ? "--am ${am_tsv}"           : ""
    def keep_arg = (keep_samples && keep_samples.size() > 0) ? "--keep ${keep_samples}" : ""
    """
    set -euo pipefail
    ${projectDir}/bin/build_carrier_matrix.py \\
        --gt ${gt} \\
        ${has_cv} ${has_ac} ${has_am} \\
        ${keep_arg} \\
        --out-long ${chunk_id}.carrier_matrix.tsv
    """
}
