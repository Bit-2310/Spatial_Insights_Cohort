# Deployment Guide for shinyapps.io

## Prerequisites
- shinyapps.io account
- `rsconnect` R package installed

## Setup

### 1. Install rsconnect
```r
install.packages("rsconnect")
```

### 2. Configure shinyapps.io account
```r
library(rsconnect)
rsconnect::setAccountInfo(
  name='<ACCOUNT>',
  token='<TOKEN>',
  secret='<SECRET>'
)
```

Get your token and secret from: https://www.shinyapps.io/admin/#/tokens

## Prepare Data for Deployment

### Important: File Size Limits
shinyapps.io has memory and file size constraints. The app is optimized to:
- Skip files larger than 30MB
- Use lazy loading for large datasets
- Cache loaded data in reactive values

### Recommended Data Structure
```
app/
├── app.R
├── R/
│   ├── load_data.R
│   ├── metrics.R
│   └── plots.R
├── data/              # Create this for deployment
│   └── processed/     # Copy SMALL sample files here
├── outputs/           # Copy summary CSVs and small plots
│   ├── qc/
│   ├── domains/
│   ├── svg/
│   └── neighborhood/
└── www/
    └── styles.css
```

### Prepare Lightweight Data
```bash
# From project root
cd app

# Create data directory
mkdir -p data/processed
mkdir -p outputs/{qc,domains,svg,neighborhood}

# Copy SMALL files only (< 30MB each)
# Check file sizes first
ls -lh ../data/processed/*.rds

# Copy small RDS files
cp ../data/processed/*.06_svg.rds data/processed/
cp ../data/processed/*.07_neighborhood.rds data/processed/

# Copy CSV summaries (always small)
cp -r ../outputs/qc outputs/
cp -r ../outputs/domains outputs/
cp -r ../outputs/svg outputs/
cp -r ../outputs/neighborhood outputs/
cp ../outputs/run_metadata.json outputs/

# Remove large RDS files if copied
find data/processed -name "*.rds" -size +30M -delete
```

## Deploy

### Option 1: Deploy from R Console
```r
library(rsconnect)
setwd("app")  # Make sure you're in the app directory

deployApp(
  appName = "spatial-insights-cohort",
  appTitle = "Spatial Insights Cohort",
  appFiles = c(
    "app.R",
    "R/",
    "data/",
    "outputs/",
    "www/"
  ),
  launch.browser = TRUE
)
```

### Option 2: Deploy from RStudio
1. Open `app/app.R` in RStudio
2. Click the "Publish" button (blue icon in top right)
3. Select "Publish Application"
4. Choose shinyapps.io
5. Select files to deploy
6. Click "Publish"

## Troubleshooting

### Error: "Application deployment failed"
- Check file sizes: `du -sh app/*`
- Ensure total app size < 1GB
- Remove large RDS files

### Error: "Memory limit exceeded"
- Reduce max_size_mb in load_data.R (currently 30MB)
- Process fewer samples
- Use CSV summaries instead of RDS files where possible

### App loads but shows "Data not available"
- Check that data/ and outputs/ directories exist in deployed app
- Verify file paths are relative (not absolute)
- Check shinyapps.io logs for errors

## Viewing Logs
```r
library(rsconnect)
showLogs(appName = "spatial-insights-cohort", streaming = TRUE)
```

## Update Deployment
```r
library(rsconnect)
setwd("app")
deployApp(appName = "spatial-insights-cohort")
```

## Performance Tips

1. **Use CSV summaries**: Prefer CSV files over RDS for shinyapps.io
2. **Limit samples**: Deploy with 2-3 samples max for demo
3. **Pre-generate plots**: Save plots as PNG and display them instead of generating on-the-fly
4. **Disable heavy features**: Comment out features that load large RDS files

## Example: Minimal Deployment
For a lightweight demo deployment:
```bash
# Keep only CSV files and small RDS
cd app
rm -rf data/processed/*.0[1-5]_*.rds  # Remove large intermediate files
rm -rf data/processed/*_domains.rds   # Remove large domain files
rm -rf data/processed/*_norm.rds      # Remove large normalized files

# Keep only:
# - *.06_svg.rds (small, ~7KB)
# - *.07_neighborhood.rds (small, ~442B)
# - All CSV files
# - PNG plots
```

## Cost Considerations
- **Free tier**: 25 active hours/month, 5 apps
- **Starter plan**: $9/month, 500 hours, better performance
- Monitor usage at: https://www.shinyapps.io/admin/#/dashboard

## Alternative: Docker Deployment
For larger datasets, consider deploying with Docker on:
- AWS ECS
- Google Cloud Run
- DigitalOcean App Platform
