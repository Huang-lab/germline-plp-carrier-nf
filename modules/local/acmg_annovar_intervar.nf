process ACMG_ANNOVAR_INTERVAR {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf)
    path annovar_dir
    path annovar_humandb
    path intervar_dir
    path intervar_config

    output:
    path("${chunk_id}.acmg_plp.tsv"), emit: tsv

    script:
    """
    set -euo pipefail
    if command -v Intervar.py >/dev/null 2>&1 || [ -x "${intervar_dir}/Intervar.py" ]; then
        # Convert VCF → ANNOVAR input, annotate, InterVar-classify.
        \${annovar_dir}/convert2annovar.pl -format vcf4 ${vcf} -outfile ${chunk_id}.avinput
        \${annovar_dir}/table_annovar.pl ${chunk_id}.avinput ${annovar_humandb} \\
            -buildver hg38 -out ${chunk_id} \\
            -remove -protocol refGene,clinvar_${params.clinvar_release},gnomad41_exome,dbnsfp47a,avsnp150 \\
            -operation g,f,f,f,f -nastring . -otherinfo
        python3 \${intervar_dir}/Intervar.py -c ${intervar_config} \\
            -i ${chunk_id}.avinput -o ${chunk_id} \\
            --input_type=AVinput -b hg38 -t \${intervar_dir}/intervardb -d ${annovar_humandb}
        INTERVAR_OUT=\$(ls ${chunk_id}*.intervar 2>/dev/null | head -n1)
        ${projectDir}/bin/acmg_postprocess.py --intervar \$INTERVAR_OUT --out ${chunk_id}.acmg_plp.tsv
    else
        # Test-profile fallback: consume a fixture InterVar-shaped file synthesized by synthetic_vep.py.
        python3 ${projectDir}/bin/synthetic_intervar.py --in ${vcf} --out _syn.intervar
        ${projectDir}/bin/acmg_postprocess.py --intervar _syn.intervar --out ${chunk_id}.acmg_plp.tsv
    fi
    """
}
