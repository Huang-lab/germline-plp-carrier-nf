#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NORM_QC }                from './modules/local/norm_qc.nf'
include { VEP_ANNOTATE }           from './modules/local/vep_annotate.nf'
include { VALIDATE_CHUNK }         from './modules/local/validate_chunk.nf'
include { CLINVAR_CLASSIFY }       from './modules/local/clinvar_classify.nf'
include { ALPHAMISSENSE_CLASSIFY } from './modules/local/alphamissense_classify.nf'
include { ACMG_ANNOVAR_INTERVAR }  from './modules/local/acmg_annovar_intervar.nf'
include { PER_GENE_QC }            from './modules/local/per_gene_qc.nf'
include { CARRIER_MATRIX }         from './modules/local/carrier_matrix.nf'
include { MANIFEST }               from './modules/local/manifest.nf'


workflow {
    if (!params.input_vcfs)      { error "Missing required param: --input_vcfs" }
    if (!params.reference_fasta) { error "Missing required param: --reference_fasta" }

    def optional_path = { p -> p ? file(p) : file("${projectDir}/tests/synthetic/EMPTY_PLACEHOLDER") }

    Channel.fromPath(params.input_vcfs, checkIfExists: true)
        .map { f -> tuple(f.baseName.replaceAll(/\.vcf(\.gz)?$/, ''), f) }
        .set { chunk_ch }

    NORM_QC(chunk_ch, file(params.reference_fasta))

    VEP_ANNOTATE(
        NORM_QC.out.vcf,
        optional_path(params.vep_cache_dir),
        file(params.reference_fasta),
        optional_path(params.clinvar_vcf),
        optional_path(params.gnomad_vcf),
        optional_path(params.am_data_tsv),
        optional_path(params.am_plugin_pm),
        optional_path(params.loftee_data_dir),
    )

    VALIDATE_CHUNK(VEP_ANNOTATE.out.vcf)

    CLINVAR_CLASSIFY(VALIDATE_CHUNK.out.vcf)
    ALPHAMISSENSE_CLASSIFY(VALIDATE_CHUNK.out.vcf, optional_path(params.am_calibration_tsv))
    ACMG_ANNOVAR_INTERVAR(
        VALIDATE_CHUNK.out.vcf,
        optional_path(params.annovar_dir),
        optional_path(params.annovar_humandb),
        optional_path(params.intervar_dir),
        optional_path(params.intervar_config),
    )

    // Per-chunk join keeps each chunk's classifications together.
    def per_chunk = VALIDATE_CHUNK.out.vcf
        .join(CLINVAR_CLASSIFY.out.tsv.map { f -> tuple(f.baseName.replaceAll(/\.clinvar_plp$/, ''), f) })
        .join(ACMG_ANNOVAR_INTERVAR.out.tsv.map { f -> tuple(f.baseName.replaceAll(/\.acmg_plp$/, ''), f) })
        .join(ALPHAMISSENSE_CLASSIFY.out.tsv.map { f -> tuple(f.baseName.replaceAll(/\.am_plp$/, ''), f) })

    // Split into typed channels for downstream processes.
    per_chunk.map { id, vcf, cv, ac, am -> tuple(id, vcf) }.set { pc_vcf }
    per_chunk.map { id, vcf, cv, ac, am -> cv }.set { pc_cv }
    per_chunk.map { id, vcf, cv, ac, am -> ac }.set { pc_ac }
    per_chunk.map { id, vcf, cv, ac, am -> am }.set { pc_am }

    PER_GENE_QC(pc_vcf, pc_cv, pc_ac, pc_am)

    def keep = params.keep_samples ? file(params.keep_samples) : file("${projectDir}/tests/synthetic/EMPTY_PLACEHOLDER")
    CARRIER_MATRIX(pc_vcf, pc_cv, pc_ac, pc_am, keep)

    MANIFEST(
        workflow.commitId ?: '',
        params.vep_cache_version ?: '',
        params.clinvar_release ?: '',
        params.gnomad_version ?: '',
    )
}
