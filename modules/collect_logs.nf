process COLLECT_LOGS {
  tag "collect_logs"
  input:
    val trigger
  output:
    path("step_summaries_long.csv"), emit: long
    path("step_summaries_wide.csv"), emit: wide
  script:
  """
  PIPELINE_ROOT=${projectDir} Rscript ${projectDir}/scripts/helpers/collect_step_summaries.R \
    --outdir ${params.outdir}
  mkdir -p ${params.outdir}/outputs/logs
  cp ${params.outdir}/outputs/logs/step_summaries_long.csv step_summaries_long.csv
  cp ${params.outdir}/outputs/logs/step_summaries_wide.csv step_summaries_wide.csv
  """
}
