process ACMG_ANNOVAR_INTERVAR {
    tag "$chunk_id"

    input:
    tuple val(chunk_id), path(vcf)
    path annovar_dir
    path annovar_humandb
    path intervar_dir
    path intervar_config

    output:
    path("${chunk_id}.acmg_plp.tsv"),        emit: tsv
    path("${chunk_id}.hg38_multianno.txt"),  emit: annovar_multianno
    path("${chunk_id}*.intervar"),           emit: intervar_raw

    script:
    def allow_stub = (params.allow_stub_acmg as boolean)
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
    elif [ "${allow_stub}" = "true" ]; then
        # Synthetic stub — ONLY for the test profile (params.allow_stub_acmg=true).
        # Produces the declared output files so DAG wiring/publishing is exercised.
        python3 ${projectDir}/bin/synthetic_intervar.py --in ${vcf} --out ${chunk_id}.hg38_multianno.txt.intervar
        printf '#Chr\\tStart\\tEnd\\tRef\\tAlt\\tFunc.refGene\\tGene.refGene\\n' > ${chunk_id}.hg38_multianno.txt
        ${projectDir}/bin/acmg_postprocess.py --intervar ${chunk_id}.hg38_multianno.txt.intervar --out ${chunk_id}.acmg_plp.tsv
    else
        echo "ERROR: ACMG requested but ANNOVAR/InterVar is not available." >&2
        echo "  Run setup/setup_annovar_intervar.sh and set params annovar_dir, annovar_humandb," >&2
        echo "  intervar_dir, intervar_config, and container_annovar_intervar." >&2
        echo "  (For pipeline testing only, set --allow_stub_acmg true to emit a synthetic stub.)" >&2
        exit 1
    fi
    """
}
