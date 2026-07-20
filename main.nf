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
include { CONCAT_CARRIER_MATRIX }  from './modules/local/concat_carrier_matrix.nf'
include { MANIFEST }               from './modules/local/manifest.nf'


// --- helpers -----------------------------------------------------------------
def parseClassifiers(spec) {
    if (spec == null || spec == '' || (spec instanceof List && spec.isEmpty()))
        return ['clinvar','acmg','am']
    def parts = (spec instanceof List) ? spec : spec.toString().split(',').collect{ it.trim().toLowerCase() }
    parts = parts.findAll{ it }
    def valid = ['clinvar','acmg','am']
    parts.each { if (!valid.contains(it)) error "Unknown classifier: '${it}'. Valid: ${valid}" }
    return valid.findAll{ parts.contains(it) }
}

// Optional-file sentinel: an empty file that build_carrier_matrix ignores.
// Avoids the "placeholder text read as IDs" bug of prior EMPTY_PLACEHOLDER.
def optionalFile(p) {
    if (p) return file(p)
    def sentinel = file("${workflow.workDir}/.empty_input")
    if (!sentinel.exists()) { sentinel.text = '' }
    return sentinel
}


workflow {
    if (!params.input_vcfs)      { error "Missing required param: --input_vcfs" }
    if (!params.reference_fasta) { error "Missing required param: --reference_fasta" }

    def classifiers = parseClassifiers(params.classifiers)
    log.info "germline-plp-carrier-nf: classifiers = ${classifiers}"

    Channel.fromPath(params.input_vcfs, checkIfExists: true)
        .map { f -> tuple(f.baseName.replaceAll(/\.vcf(\.gz)?$/, ''), f) }
        .set { chunk_ch }

    NORM_QC(chunk_ch, file(params.reference_fasta))

    VEP_ANNOTATE(
        NORM_QC.out.vcf,
        optionalFile(params.vep_cache_dir),
        file(params.reference_fasta),
        optionalFile(params.clinvar_vcf),
        optionalFile(params.gnomad_vcf),
        optionalFile(params.am_data_tsv),
        optionalFile(params.am_plugin_pm),
        optionalFile(params.loftee_data_dir),
    )

    VALIDATE_CHUNK(VEP_ANNOTATE.out.vcf, classifiers.join(','))

    // Empty channel emitting an empty file per chunk — used when a classifier is off.
    def empty_ch = { VALIDATE_CHUNK.out.vcf.map { id, vcf -> file("${workflow.workDir}/.empty_${id}") } }

    def clinvar_tsv_ch = classifiers.contains('clinvar')
        ? CLINVAR_CLASSIFY(VALIDATE_CHUNK.out.vcf).tsv
        : empty_ch()
    def am_tsv_ch = classifiers.contains('am')
        ? ALPHAMISSENSE_CLASSIFY(VALIDATE_CHUNK.out.vcf, optionalFile(params.am_calibration_tsv)).tsv
        : empty_ch()
    def acmg_tsv_ch = classifiers.contains('acmg')
        ? ACMG_ANNOVAR_INTERVAR(
              VALIDATE_CHUNK.out.vcf,
              optionalFile(params.annovar_dir),
              optionalFile(params.annovar_humandb),
              optionalFile(params.intervar_dir),
              optionalFile(params.intervar_config),
          ).tsv
        : empty_ch()

    // Per-chunk join: tuple(chunk_id, vcf, clinvar_tsv, acmg_tsv, am_tsv)
    // Each classifier TSV filename encodes the chunk_id as basename prefix.
    def keyByChunk = { ch, suffix ->
        ch.map { f -> tuple(f.baseName.replaceAll(suffix, ''), f) }
    }
    def per_chunk = VALIDATE_CHUNK.out.vcf
        .join(keyByChunk(clinvar_tsv_ch, /\.clinvar_plp$/))
        .join(keyByChunk(acmg_tsv_ch,    /\.acmg_plp$/))
        .join(keyByChunk(am_tsv_ch,      /\.am_plp$/))

    per_chunk.map { id, vcf, cv, ac, am -> tuple(id, vcf) }.set { pc_vcf }
    per_chunk.map { id, vcf, cv, ac, am -> cv }.set { pc_cv }
    per_chunk.map { id, vcf, cv, ac, am -> ac }.set { pc_ac }
    per_chunk.map { id, vcf, cv, ac, am -> am }.set { pc_am }

    PER_GENE_QC(pc_vcf, pc_cv, pc_ac, pc_am, classifiers.join(','))

    CARRIER_MATRIX(pc_vcf, pc_cv, pc_ac, pc_am,
                   optionalFile(params.keep_samples), classifiers.join(','))

    CONCAT_CARRIER_MATRIX(CARRIER_MATRIX.out.long_tsv.collect())

    MANIFEST(
        workflow.commitId ?: '',
        params.vep_cache_version ?: '',
        params.clinvar_release ?: '',
        params.gnomad_version ?: '',
        classifiers.join(','),
    )
}
