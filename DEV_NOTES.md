# Development Notes

## Environment
Use conda to create the R environment:

```
conda env create -f environment.yml
conda activate spatial_insights_cohort
```

## Manifest
The manifest must include `sample_id`. You can optionally add `spe_source` for row-level input overrides.

## Run locally with Nextflow
From the repo root:

```
nextflow run main.nf -profile local \
  --manifest assets/cohort_manifest_template.csv \
  --steps "01,02,03,04,05,06,07" \
  --max_samples 12 \
  -resume
```

## Validate only
Run validation without executing any steps:

```
nextflow run main.nf -profile local \
  --manifest assets/cohort_manifest_template.csv \
  --steps "01,02,03,04,05,06,07" \
  --validate_only true \
  --validate_strict true
```

## Using an input SpatialExperiment RDS
```
nextflow run main.nf -profile local \
  --manifest assets/cohort_manifest_template.csv \
  --dataset dlpfc \
  --input_spe_rds /path/to/cohort_spe.rds \
  --steps "01,02,03,04,05,06,07" \
  -resume
```

## Key parameters
1. `--dataset` (default: `dlpfc`)
2. `--input_spe_rds` (optional)
3. `--validate_only` (default: `false`)
4. `--validate_strict` (default: `true`)
5. `--steps` (comma-separated step ids)

## Outputs
1. Caches: `data/processed/`
2. QC: `outputs/qc/`
3. Integration: `outputs/integration/`
4. Domains: `outputs/domains/`
5. SVG: `outputs/svg/`
6. Neighborhood: `outputs/neighborhood/`
7. Validation: `outputs/validation/`
8. Run metadata: `outputs/run_metadata.json`, `outputs/run_metadata.txt`
9. Nextflow reports: `outputs/nf/`
10. Step summaries: `outputs/<step>/summary/<sample_id>.tsv`
11. Aggregated logs: `outputs/logs/step_summaries_long.csv`, `outputs/logs/step_summaries_wide.csv`
