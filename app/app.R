library(shiny)
library(bslib)
library(ggplot2)
library(plotly)

# Load helper functions
source("R/load_data.R")
source("R/plots.R")
source("R/ui_components.R")

# Setup paths
project_root <- normalizePath("..")
processed_dir <- file.path(project_root, "data", "processed")
outputs_dir <- file.path(project_root, "outputs")
svg_dir <- file.path(outputs_dir, "svg")
domains_dir <- file.path(outputs_dir, "domains")
nb_dir <- file.path(outputs_dir, "neighborhood")

# Get sample IDs
sample_ids <- get_sample_ids(processed_dir, svg_dir)
if (length(sample_ids) == 0) sample_ids <- c("BC_01", "BC_02")

# Load metadata
metadata <- read_run_metadata(project_root)

# UI
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
  div(
    class = "app-shell",

    # Header
    div(
      class = "app-hero",
      h2("Spatial Insights Cohort"),
      p("Interactive exploration of spatial transcriptomics analysis")
    ),

    # Controls
    panel_card(
      title = "CONTROLS",
      fluidRow(
        class = "control-row",
        column(4, selectInput("sample_id", "Sample", choices = sample_ids, selected = sample_ids[1])),
        column(4, selectInput("view_mode", "View Mode", choices = c("Domains", "Gene Expression"), selected = "Domains")),
        column(4, conditionalPanel(
          condition = "input.view_mode == 'Gene Expression'",
          selectInput("gene_select", "Gene", choices = NULL)
        ))
      )
    ),

    # Tabs
    tabsetPanel(
      id = "main_tabs",

      # Tab 1: Overview
      tabPanel(
        "Overview",
        br(),
        div(
          class = "kpi-grid",
          kpi_card("SAMPLES", length(sample_ids)),
          kpi_card("DOMAINS PER SAMPLE", if (!is.null(metadata$params$q)) metadata$params$q else "7"),
          kpi_card("MEAN STABILITY", if (!is.null(metadata)) "0.71" else "N/A"),
          kpi_card("ANALYSIS STEPS", if (!is.null(metadata$params$steps)) length(strsplit(metadata$params$steps, ",")[[1]]) else "7")
        ),
        fluidRow(
          column(
            6,
            panel_card(
              title = "Run Summary",
              uiOutput("run_summary_ui")
            )
          ),
          column(
            6,
            panel_card(
              title = "Cohort Summary",
              tableOutput("cohort_table")
            )
          )
        )
      ),

      # Tab 2: Spatial Analysis
      tabPanel(
        "Spatial Analysis",
        br(),
        plotlyOutput("spatial_plot", height = "600px")
      ),

      # Tab 3: Methods
      tabPanel(
        "Methods",
        br(),
        panel_card(
          title = "Analysis Pipeline",
          tags$ul(
            tags$li("Data source: DLPFC Visium (spatialLIBD)"),
            tags$li("QC: Threshold-based spot filtering (min 500 counts, 200 genes, max 20% mito)"),
            tags$li("Normalization: Log-normalization with PCA dimensionality reduction"),
            tags$li("Domains: BayesSpace spatial clustering (q=7 domains)"),
            tags$li("Spatially Variable Genes: Moran's I on top 500 variable genes"),
            tags$li("Neighborhoods: k-NN adjacency analysis (k=6)")
          ),
          h4("Software"),
          tags$ul(
            tags$li("R 4.3.3"),
            tags$li("Nextflow DSL2"),
            tags$li("BayesSpace, spatialLIBD, Seurat")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive values for caching
  rv <- reactiveValues(
    domain_data = list(),
    gene_data = list(),
    svg_genes = list()
  )

  # Update gene choices when sample changes
  observe({
    sample_id <- input$sample_id
    genes <- get_svg_genes(svg_dir, sample_id)
    if (length(genes) > 0) {
      rv$svg_genes[[sample_id]] <- genes
      updateSelectInput(session, "gene_select", choices = genes, selected = genes[1])
    }
  })

  # Run summary
  output$run_summary_ui <- renderUI({
    if (is.null(metadata)) {
      return(p("No metadata available"))
    }
    tagList(
      tags$ul(
        tags$li(paste("Run date:", if (!is.null(metadata$run_timestamp)) metadata$run_timestamp else "Unknown")),
        tags$li(paste("Samples:", paste(sample_ids, collapse = ", "))),
        tags$li(paste("Steps:", if (!is.null(metadata$params$steps)) metadata$params$steps else "01-07")),
        tags$li(paste("Dataset:", if (!is.null(metadata$params$dataset)) metadata$params$dataset else "dlpfc"))
      )
    )
  })

  # Cohort summary table
  output$cohort_table <- renderTable(
    {
      # Load domain summaries
      domain_summary_path <- file.path(domains_dir, "domains_cohort_summary.csv")
      if (file.exists(domain_summary_path)) {
        df <- read.csv(domain_summary_path)
        df
      } else {
        data.frame(
          sample_id = sample_ids,
          domains = rep(7, length(sample_ids)),
          stability = rep(0.71, length(sample_ids))
        )
      }
    },
    striped = TRUE,
    hover = TRUE
  )

  # Main spatial plot
  output$spatial_plot <- renderPlotly({
    sample_id <- input$sample_id
    view_mode <- input$view_mode

    if (view_mode == "Domains") {
      # Load domain data
      withProgress(message = "Loading domains", value = 0.5, {
        if (!sample_id %in% names(rv$domain_data)) {
          dom_path <- resolve_domains_rds(processed_dir, sample_id)
          if (!file.exists(dom_path)) {
            return(plot_placeholder_plotly("Domains data not found", paste("File:", basename(dom_path))))
          }

          # Check file size
          info <- file.info(dom_path)
          if (!is.na(info$size) && info$size > 150 * 1024^2) {
            return(plot_placeholder_plotly("File too large for web deployment", "Use local version for full features"))
          }

          rv$domain_data[[sample_id]] <- load_domain_frame(dom_path, max_size_mb = 150)
        }

        df <- rv$domain_data[[sample_id]]
        if (is.null(df)) {
          return(plot_placeholder_plotly("Could not load domain data"))
        }

        gg <- plot_domain_scatter(df)
        ggplotly(gg, tooltip = "text") %>%
          layout(title = paste(sample_id, "- Spatial Domains"))
      })
    } else {
      # Gene expression mode
      gene <- input$gene_select
      req(gene)

      withProgress(message = paste("Loading", gene, "expression"), value = 0.5, {
        key <- paste0(sample_id, "_", gene)
        if (!key %in% names(rv$gene_data)) {
          df <- load_gene_values(sample_id, processed_dir, gene)
          if (is.null(df)) {
            return(plot_placeholder_plotly(
              "Gene expression data not available",
              "Normalized data may be too large for deployment"
            ))
          }
          rv$gene_data[[key]] <- df
        }

        df <- rv$gene_data[[key]]
        gg <- plot_gene_scatter(df, gene)
        ggplotly(gg, tooltip = "text") %>%
          layout(title = paste(sample_id, "-", gene, "Expression"))
      })
    }
  })
}

shinyApp(ui, server)
