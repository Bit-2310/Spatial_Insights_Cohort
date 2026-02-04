process INTEGRATE {
  tag "cohort"
  input:
    path(norm_files)
  output:
    path("cohort.04_integrate.rds")
  script:
  """
  PIPELINE_ROOT=${projectDir} Rscript ${projectDir}/scripts/04_integrate.R \
    --in_files ${norm_files.join(',')} \
    --out_processed ${params.outdir}/data/processed \
    --out_outputs ${params.outdir}/outputs \
    --log_file ${params.outdir}/logs/04_integrate/cohort.log
  ln -sf ${params.outdir}/data/processed/cohort.04_integrate.rds cohort.04_integrate.rds
  """
}
