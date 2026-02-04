process INGEST {
  tag { sample_id }
  input:
    tuple val(sample_id), val(spe_source)
  output:
    tuple val(sample_id), path("${sample_id}.01_ingest.rds")
  script:
  def spe_arg = spe_source ? "--spe_source ${spe_source}" : ""
  def input_arg = params.input_spe_rds ? "--input_spe_rds ${params.input_spe_rds}" : ""
  """
  PIPELINE_ROOT=${projectDir} Rscript ${projectDir}/scripts/01_ingest.R \
    --sample_id ${sample_id} \
    --dataset ${params.dataset} \
    ${input_arg} \
    ${spe_arg} \
    --out_processed ${params.outdir}/data/processed \
    --out_outputs ${params.outdir}/outputs \
    --log_file ${params.outdir}/logs/01_ingest/${sample_id}.log
  ln -sf ${params.outdir}/data/processed/${sample_id}.01_ingest.rds ${sample_id}.01_ingest.rds
  """
}
