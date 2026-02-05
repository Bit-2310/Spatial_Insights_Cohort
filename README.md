# Spatial Insights Cohort

This is a friendly, end-to-end spatial transcriptomics pipeline in R to process a cohort of human DLPFC Visium samples, produce consistent cohort summaries, and power a Shiny dashboard.

## What you get
After a run, you will have clean per-sample outputs plus cohort-level summaries that make it easy to validate data quality, compare samples, and communicate results with stakeholders.

**Core results**
- QC summaries and plots to confirm data health and filtering impact
- Normalized assays with PCA diagnostics for each sample
- Spatial domains with stability metrics and per-sample domain counts
- Spatially variable genes (SVGs) with ranked gene lists and cohort consensus
- Neighborhood adjacency summaries to compare spatial organization across samples

**Business-ready artifacts**
- CSV summaries for downstream analytics or reporting
- Plots for quick inspection or slide decks
- A Shiny dashboard for interactive review

## Data
Right now the pipeline uses the DLPFC Visium cohort from spatialLIBD by default. In the future, this will be expanded to support additional datasets and platforms.

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

## Pipeline flow (high level)
1. Ingest and cache each sample
2. QC filtering with clear thresholds
3. Normalization and per-sample PCA
4. Spatial domains with stability checks
5. Spatially variable genes
6. Neighborhood adjacency summaries

## Outputs
Outputs are organized under `outputs/` and cached data under `data/processed/`.

**Key folders**
- `outputs/qc`: QC summaries, thresholds, and plots
- `outputs/integration`: Normalization summaries and (optional) cohort integration artifacts
- `outputs/domains`: Domain summaries, stability, and plots
- `outputs/svg`: Ranked SVGs and cohort consensus lists
- `outputs/neighborhood`: Neighborhood matrices and cohort summaries
- `outputs/validation`: Run validation reports
- `outputs/run_metadata.json`: Run metadata and provenance

## Project status
Project status: stable. No new features planned.
