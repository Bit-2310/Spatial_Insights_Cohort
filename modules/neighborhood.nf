process NEIGHBORHOOD {
  tag { sample_id }
  input:
    tuple val(sample_id), path(in_rds)
  output:
    tuple val(sample_id), path("${sample_id}.07_neighborhood.rds")
  script:
  """
  PIPELINE_ROOT=${projectDir} Rscript ${projectDir}/scripts/07_neighborhood.R \
    --in_rds ${in_rds} \
    --out_processed ${params.outdir}/data/processed \
    --out_outputs ${params.outdir}/outputs \
    --log_file ${params.outdir}/logs/07_neighborhood/${sample_id}.log
  ln -sf ${params.outdir}/data/processed/${sample_id}.07_neighborhood.rds ${sample_id}.07_neighborhood.rds
  """
}
