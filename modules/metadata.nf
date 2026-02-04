process WRITE_METADATA {
  tag "run_metadata"
  input:
    val sample_ids
  output:
    path("run_metadata.json"), emit: json
    path("run_metadata.txt"), emit: txt
  script:
  def ids = (sample_ids instanceof List) ? sample_ids.join(',') : sample_ids
  def input_arg = params.input_spe_rds ? "--input_spe_rds ${params.input_spe_rds}" : ""
  """
  PIPELINE_ROOT=${projectDir} Rscript ${projectDir}/scripts/helpers/write_run_metadata.R \
    --out_json run_metadata.json \
    --out_txt run_metadata.txt \
    --manifest ${params.manifest} \
    --max_samples ${params.max_samples} \
    --steps ${params.steps} \
    --outdir ${params.outdir} \
    --dataset ${params.dataset} \
    ${input_arg} \
    --qc_min_counts 500 \
    --qc_min_genes 200 \
    --qc_max_mito 20 \
    --qc_keep_tissue_only TRUE \
    --q_sweep 6,7,8 \
    --hvg_count_norm 2000 \
    --hvg_count_svg 500 \
    --knn_k_svg 6 \
    --knn_k_nb 7 \
    --sample_ids ${ids}
  mkdir -p ${params.outdir}/outputs
  cp run_metadata.json ${params.outdir}/outputs/run_metadata.json
  cp run_metadata.txt ${params.outdir}/outputs/run_metadata.txt
  """
}
