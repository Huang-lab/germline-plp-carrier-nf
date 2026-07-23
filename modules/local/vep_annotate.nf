process VEP_ANNOTATE {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf)
    path vep_cache_dir
    path reference_fasta
    path clinvar_vcf
    path gnomad_vcf
    path am_data_tsv
    path am_plugin_pm
    path loftee_data_dir

    output:
    tuple val(chunk_id), path("${chunk_id}.vep.vcf.gz"), emit: vcf

    script:
    def cache_ver = params.vep_cache_version ?: params.vep_assembly
    """
    set -euo pipefail
    if command -v vep >/dev/null 2>&1; then
        # --custom VCF tracks need a tabix index adjacent to the file. Nextflow
        # stages the .vcf.gz but not its index, so build one in the task dir if
        # absent (tabix ships with the VEP container's htslib).
        for f in "${clinvar_vcf}" "${gnomad_vcf}"; do
            if [ -s "\$f" ] && [ ! -e "\$f.tbi" ] && [ ! -e "\$f.csi" ]; then
                tabix -p vcf "\$f" || true
            fi
        done

        # Build optional tracks/plugins ONLY for references actually supplied.
        # Unused optional inputs are 0-byte placeholders (see main.nf optionalFile),
        # so a size test cleanly enables/disables each one.
        EXTRA=""
        if [ -s "${clinvar_vcf}" ]; then
            EXTRA="\$EXTRA --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNSIGCONF"
        fi
        if [ -s "${gnomad_vcf}" ]; then
            EXTRA="\$EXTRA --custom ${gnomad_vcf},gnomAD,vcf,exact,0,AF,${params.gnomad_popmax_field}"
        fi
        if [ -s "${am_data_tsv}" ]; then
            EXTRA="\$EXTRA --dir_plugins \$(dirname \$(readlink -f ${am_plugin_pm})) --plugin AlphaMissense,file=${am_data_tsv}"
        fi
        if [ -s "${loftee_data_dir}/human_ancestor.fa.gz" ]; then
            EXTRA="\$EXTRA --plugin LoF,loftee_path:/opt/loftee,human_ancestor_fa:${loftee_data_dir}/human_ancestor.fa.gz"
        fi

        vep --input_file ${vcf} --output_file ${chunk_id}.vep.vcf.gz --vcf \\
            --compress_output bgzip \\
            --force_overwrite --no_stats --offline --cache \\
            --dir_cache ${vep_cache_dir} --cache_version ${cache_ver} \\
            --species ${params.vep_species} --assembly ${params.vep_assembly} \\
            --fasta ${reference_fasta} \\
            --symbol --canonical --mane \\
            --pick --pick_order mane_select,canonical,tsl,biotype,rank,ccds,length \\
            \$EXTRA
    else
        # Test-profile fallback: synthesize a minimal VEP-like CSQ header + a CSQ INFO tag
        # for each variant, then bgzip-equivalent gzip.
        python3 ${projectDir}/bin/synthetic_vep.py --in ${vcf} --out ${chunk_id}.vep.vcf
        gzip -f ${chunk_id}.vep.vcf
    fi
    """
}
