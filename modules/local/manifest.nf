process MANIFEST {
    input:
    val pipeline_sha
    val vep_cache_version
    val clinvar_release
    val gnomad_version

    output:
    path("manifest.json")

    script:
    """
    ${projectDir}/bin/make_manifest.py \\
        --pipeline-sha "${pipeline_sha ?: ''}" \\
        --reference-build ${params.vep_assembly} \\
        --vep-cache "${vep_cache_version ?: ''}" \\
        --clinvar-release "${clinvar_release ?: ''}" \\
        --gnomad-version "${gnomad_version ?: ''}" \\
        --am-calibration-version "user-supplied" \\
        --out manifest.json
    """
}
