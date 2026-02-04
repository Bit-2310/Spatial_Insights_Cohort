process VALIDATE {
  tag "validate"
  input:
    val sample_ids
    val trigger
  output:
    path("validation_report.json"), emit: json
    path("validation_report.txt"), emit: txt
    val(true), emit: done
  script:
  def ids = (sample_ids instanceof List) ? sample_ids.join(',') : sample_ids
  """
  PIPELINE_ROOT=${projectDir} Rscript ${projectDir}/scripts/helpers/validate_outputs.R \
    --outdir ${params.outdir} \
    --manifest ${params.manifest} \
    --steps ${params.steps} \
    --sample_ids ${ids} \
    --strict ${params.validate_strict}
  mkdir -p ${params.outdir}/outputs/validation
  cp ${params.outdir}/outputs/validation/validation_report.json validation_report.json
  cp ${params.outdir}/outputs/validation/validation_report.txt validation_report.txt
  """
}
