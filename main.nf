#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NORM_QC }                from './modules/local/norm_qc.nf'
include { VEP_ANNOTATE }           from './modules/local/vep_annotate.nf'
include { VALIDATE_CHUNK }         from './modules/local/validate_chunk.nf'
include { CLINVAR_CLASSIFY }       from './modules/local/clinvar_classify.nf'
include { ALPHAMISSENSE_CLASSIFY } from './modules/local/alphamissense_classify.nf'
include { ACMG_ANNOVAR_INTERVAR }  from './modules/local/acmg_annovar_intervar.nf'
include { PER_GENE_QC }            from './modules/local/per_gene_qc.nf'
include { CARRIER_GT }             from './modules/local/carrier_gt.nf'
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

// Optional-file sentinel: a UNIQUELY-NAMED empty file per logical input, so a
// process taking several optional paths never sees two files with the same
// name (Nextflow rejects staging duplicate basenames into one task dir).
def optionalFile(p, tag) {
    if (p) return file(p)
    def sentinel = file("${workflow.workDir}/.empty__${tag}")
    if (!sentinel.exists()) { sentinel.text = '' }
    return sentinel
}


// A required path param must be a non-empty string. A Boolean (e.g. `--flag`
// with no value) or null is an error with a clear message.
def requirePath(name, val) {
    if (val == null || val instanceof Boolean || val.toString().trim() == '') {
        error "Param --${name} must be a path/glob string, got: ${val}. " +
              "If launching via LSF, the value may have been lost — pass the " +
              "literal path (do not rely on an un-exported shell variable)."
    }
    return val.toString()
}

// QC-only entry point: run with `-entry NORM_ONLY` to produce and publish
// norm_qc/*.norm.vcf.gz and stop. No VEP, no classifiers.
workflow NORM_ONLY {
    requirePath('input_vcfs', params.input_vcfs)
    requirePath('reference_fasta', params.reference_fasta)

    Channel.fromPath(params.input_vcfs, checkIfExists: true)
        .map { f -> tuple(f.baseName.replaceAll(/\.vcf(\.gz)?$/, ''), f) }
        .set { chunk_ch }

    NORM_QC(chunk_ch, file(params.reference_fasta))
}

workflow {
    requirePath('input_vcfs', params.input_vcfs)
    requirePath('reference_fasta', params.reference_fasta)

    def classifiers = parseClassifiers(params.classifiers)
    log.info "germline-plp-carrier-nf: classifiers = ${classifiers}"

    Channel.fromPath(params.input_vcfs, checkIfExists: true)
        .map { f -> tuple(f.baseName.replaceAll(/\.vcf(\.gz)?$/, ''), f) }
        .set { chunk_ch }

    NORM_QC(chunk_ch, file(params.reference_fasta))

    VEP_ANNOTATE(
        NORM_QC.out.vcf,
        optionalFile(params.vep_cache_dir, 'vep_cache_dir'),
        file(params.reference_fasta),
        optionalFile(params.clinvar_vcf, 'clinvar_vcf'),
        optionalFile(params.gnomad_vcf, 'gnomad_vcf'),
        optionalFile(params.am_data_tsv, 'am_data_tsv'),
        optionalFile(params.am_plugin_pm, 'am_plugin_pm'),
        optionalFile(params.loftee_data_dir, 'loftee_data_dir'),
    )

    VALIDATE_CHUNK(VEP_ANNOTATE.out.vcf, classifiers.join(','))

    // When a classifier is OFF, emit a per-chunk EMPTY placeholder named exactly
    // like that classifier's real output (<id>.<suffix>), so (a) it keys to the
    // chunk in the join below and (b) its basename never collides with another
    // classifier's placeholder. Downstream steps skip zero-byte tables.
    def emptyPerChunk = { suffix -> VALIDATE_CHUNK.out.vcf.map { id, vcf ->
        def f = file("${workflow.workDir}/${id}.${suffix}")
        if (!f.exists()) { f.text = '' }
        return f
    } }

    def clinvar_tsv_ch = classifiers.contains('clinvar')
        ? CLINVAR_CLASSIFY(VALIDATE_CHUNK.out.vcf).tsv
        : emptyPerChunk('clinvar_plp')
    def am_tsv_ch = classifiers.contains('am')
        ? ALPHAMISSENSE_CLASSIFY(VALIDATE_CHUNK.out.vcf, optionalFile(params.am_calibration_tsv, 'am_calibration_tsv')).tsv
        : emptyPerChunk('am_plp')
    def acmg_tsv_ch = classifiers.contains('acmg')
        ? ACMG_ANNOVAR_INTERVAR(
              VALIDATE_CHUNK.out.vcf,
              optionalFile(params.annovar_dir, 'annovar_dir'),
              optionalFile(params.annovar_humandb, 'annovar_humandb'),
              optionalFile(params.intervar_dir, 'intervar_dir'),
              optionalFile(params.intervar_config, 'intervar_config'),
          ).tsv
        : emptyPerChunk('acmg_plp')

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

    // Efficient GT extraction (bcftools, C-speed) then pure-Python assembly.
    CARRIER_GT(per_chunk)   // input tuple(id, vcf, cv, ac, am) → tuple(id, gt)

    def cm_in = per_chunk
        .map { id, vcf, cv, ac, am -> tuple(id, cv, ac, am) }
        .join(CARRIER_GT.out.gt)   // tuple(id, cv, ac, am, gt)

    CARRIER_MATRIX(cm_in,
                   optionalFile(params.keep_samples, 'keep_samples'), classifiers.join(','))

    CONCAT_CARRIER_MATRIX(CARRIER_MATRIX.out.long_tsv.collect())

    MANIFEST(
        workflow.commitId ?: '',
        params.vep_cache_version ?: '',
        params.clinvar_release ?: '',
        params.gnomad_version ?: '',
        classifiers.join(','),
    )
}
