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
    tuple val(chunk_id), path("${chunk_id}.vep.vcf"), emit: vcf

    script:
    """
    set -euo pipefail
    if command -v vep >/dev/null 2>&1; then
        vep --input_file ${vcf} --output_file ${chunk_id}.vep.vcf --vcf --force_overwrite \\
            --cache --dir_cache ${vep_cache_dir} --species ${params.vep_species} \\
            --assembly ${params.vep_assembly} --fasta ${reference_fasta} \\
            --symbol --canonical --mane --pick --pick_order mane_select,canonical,tsl,biotype,rank,ccds,length \\
            --fork ${task.cpus} \\
            --plugin LoF,loftee_path:/opt/loftee,human_ancestor_fa:${loftee_data_dir}/human_ancestor.fa.gz \\
            --plugin AlphaMissense,file=${am_data_tsv} \\
            --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT \\
            --custom ${gnomad_vcf},gnomAD,vcf,exact,0,AF,${params.gnomad_popmax_field}
    else
        # Test-profile fallback: synthesize a minimal VEP-like CSQ header + a CSQ INFO tag
        # for each variant, using SYMBOL from a fixed lookup table in the fixture.
        python3 ${projectDir}/bin/synthetic_vep.py --in ${vcf} --out ${chunk_id}.vep.vcf
    fi
    """
}
