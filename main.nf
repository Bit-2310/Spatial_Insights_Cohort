nextflow.enable.dsl=2

params.manifest = params.manifest ?: 'assets/cohort_manifest_template.csv'
params.max_samples = params.max_samples ?: 0
params.steps = params.steps ?: '01,02,03,05,06,07'
params.outdir = params.outdir ?: projectDir
params.dataset = params.dataset ?: 'dlpfc'
params.input_spe_rds = params.input_spe_rds ?: ''
params.validate_strict = params.validate_strict == null ? true : params.validate_strict
params.validate_only = params.validate_only == null ? false : params.validate_only

// ensure output folders exist before workflow starts
new File("${params.outdir}/outputs/nf").mkdirs()
new File("${params.outdir}/logs").mkdirs()

include { INGEST } from './modules/ingest.nf'
include { QC } from './modules/qc.nf'
include { NORMALIZE } from './modules/normalize.nf'
include { INTEGRATE } from './modules/integrate.nf'
include { DOMAINS } from './modules/domains.nf'
include { SVG } from './modules/svg.nf'
include { NEIGHBORHOOD } from './modules/neighborhood.nf'
include { WRITE_METADATA } from './modules/metadata.nf'
include { VALIDATE } from './modules/validate.nf'
include { COLLECT_LOGS } from './modules/collect_logs.nf'

workflow {
  def steps = params.steps.split(',').collect { it.trim() }

  Channel
    .fromPath(params.manifest)
    .splitCsv(header:true)
    .map { row ->
      def sidKey = row.keySet().find { it.replaceAll('\"','') == 'sample_id' }
      def srcKey = row.keySet().find { it.replaceAll('\"','') == 'spe_source' }
      def sid = sidKey ? row[sidKey] : null
      def src = srcKey ? row[srcKey] : ''
      if (sid instanceof String) sid = sid.replaceAll('^\"|\"$', '')
      if (src instanceof String) src = src.replaceAll('^\"|\"$', '')
      tuple(sid, src)
    }
    .filter { sid, src -> sid != null && sid != '' }
    .ifEmpty { error "No sample_id values found in manifest" }
    .set { sample_rows_all }

  sample_rows = sample_rows_all
  if (params.max_samples && params.max_samples.toInteger() > 0) {
    sample_rows = sample_rows_all.take(params.max_samples.toInteger())
  }

  sample_ids = sample_rows.map { sid, src -> sid }
  sample_ids_list = sample_ids.collect()
  WRITE_METADATA(sample_ids_list)

  if (!params.validate_only) {
    if (steps.contains('01')) {
      ingest_out = INGEST(sample_rows)
    } else {
      ingest_out = sample_ids.map { sid -> tuple(sid, file("${params.outdir}/data/processed/${sid}.01_ingest.rds")) }
    }

    if (steps.contains('02')) {
      qc_out = QC(ingest_out)
    } else {
      qc_out = ingest_out.map { sid, path -> tuple(sid, file("${params.outdir}/data/processed/${sid}.02_qc.rds")) }
    }

    if (steps.contains('03')) {
      norm_out = NORMALIZE(qc_out)
    } else {
      norm_out = qc_out.map { sid, path -> tuple(sid, file("${params.outdir}/data/processed/${sid}.03_normalize.rds")) }
    }

    if (steps.contains('04')) {
      norm_files = norm_out.map { sid, path -> path }.collect()
      integrate_out = INTEGRATE(norm_files)
    }

    if (steps.contains('05')) {
      dom_out = DOMAINS(norm_out)
    } else {
      dom_out = norm_out.map { sid, path -> tuple(sid, file("${params.outdir}/data/processed/${sid}.05_domains.rds")) }
    }

    if (steps.contains('06')) {
      svg_out = SVG(dom_out)
    } else {
      svg_out = dom_out.map { sid, path -> tuple(sid, file("${params.outdir}/data/processed/${sid}.06_svg.rds")) }
    }

    if (steps.contains('07')) {
      nb_out = NEIGHBORHOOD(dom_out)
      if (steps.contains('06')) {
        final_done = svg_out.map { true }.mix(nb_out.map { true }).collect()
      } else {
        final_done = nb_out.map { true }.collect()
      }
    } else if (steps.contains('06')) {
      final_done = svg_out.map { true }.collect()
    } else if (steps.contains('05')) {
      final_done = dom_out.map { true }.collect()
    } else if (steps.contains('04')) {
      final_done = integrate_out.map { true }.collect()
    } else if (steps.contains('03')) {
      final_done = norm_out.map { true }.collect()
    } else if (steps.contains('02')) {
      final_done = qc_out.map { true }.collect()
    } else {
      final_done = ingest_out.map { true }.collect()
    }
  } else {
    final_done = Channel.value(true)
  }

  VALIDATE(sample_ids_list, final_done)
  COLLECT_LOGS(VALIDATE.out.done)
}
