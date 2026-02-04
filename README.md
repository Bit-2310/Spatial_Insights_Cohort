# Spatial Insights Cohort

This project runs a spatial transcriptomics pipeline in R to process a cohort of human DLPFC Visium samples and generate consistent cohort summaries plus a Shiny dashboard.

## Data
The pipeline uses the DLPFC Visium cohort from spatialLIBD by default.

## Run locally
From the repo root, run:

```
nextflow run main.nf -profile local --manifest assets/cohort_manifest_template.csv --steps "01,02,03,04,05,06,07" --max_samples 12 -resume
```

## Shiny dashboard
Run after the pipeline completes:

```
R -e "shiny::runApp('app')"
```

## Project status
Project status: stable. No new features planned.
