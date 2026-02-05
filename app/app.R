#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(ggplot2)
  library(DT)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
})

source(file.path("R", "load_data.R"))
source(file.path("R", "plots.R"))
source(file.path("R", "metrics.R"))

find_project_root <- function(start_dir, max_up = 3) {
  dir <- normalizePath(start_dir)
  for (i in seq_len(max_up + 1)) {
    if (dir.exists(file.path(dir, "data", "processed")) && dir.exists(file.path(dir, "outputs"))) {
      return(dir)
    }
    dir <- normalizePath(file.path(dir, ".."))
  }
  normalizePath(start_dir)
}

app_dir <- normalizePath(getwd())
project_root <- find_project_root(app_dir, max_up = 3)
processed_dir <- file.path(project_root, "data", "processed")
qc_path <- file.path(project_root, "outputs", "qc", "qc_summary_postfilter.csv")
integration_dir <- file.path(project_root, "outputs", "integration")
domains_dir <- file.path(project_root, "outputs", "domains")
svg_dir <- file.path(project_root, "outputs", "svg")
nb_dir <- file.path(project_root, "outputs", "neighborhood")

set_metrics_context(list(
  project_root = project_root,
  processed_dir = processed_dir,
  qc_path = qc_path,
  integration_dir = integration_dir,
  domains_dir = domains_dir,
  svg_dir = svg_dir,
  nb_dir = nb_dir
))

sample_ids <- get_sample_ids(processed_dir, svg_dir)

# Startup report
cat("Startup report:\n")
cat("- app_dir:", app_dir, "\n")
cat("- project_root:", project_root, "\n")
cat("- samples_detected:", length(sample_ids), "\n")
cat("- qc_summary:", file.exists(qc_path), "\n")
cat("- domains_summary:", file.exists(file.path(domains_dir, "domains_cohort_summary.csv")), "\n")
cat("- svg_summary:", file.exists(file.path(svg_dir, "svg_cohort_summary.csv")), "\n")
cat("- neighborhood_summary:", file.exists(file.path(nb_dir, "neighborhood_cohort_summary.csv")), "\n")
cat("- cohort_umap_plot:", file.exists(file.path(integration_dir, "plots", "cohort_umap_by_sample.png")), "\n")

  if (length(sample_ids) == 0) {
    ui <- fluidPage(
      theme = bs_theme(version = 5, base_font = font_google("IBM Plex Sans")),
      h2("Spatial Insights Cohort"),
      p("No outputs detected yet. Run the pipeline to populate the dashboard."),
      div(
        class = "panel-card",
        h4("Run pipeline stages"),
        tags$ul(
          tags$li("Ingest: Rscript scripts/01_ingest.R"),
          tags$li("QC: Rscript scripts/02_qc.R"),
          tags$li("Normalize: Rscript scripts/03_normalize.R"),
          tags$li("Domains: Rscript scripts/05_domains.R"),
          tags$li("Spatial Genes: Rscript scripts/06_svg.R"),
          tags$li("Neighborhoods: Rscript scripts/07_neighborhood.R")
        )
      )
    )
  server <- function(input, output, session) {}
  shinyApp(ui, server)
  quit(save = "no")
}

# Defaults
first_sample <- sample_ids[1]
second_sample <- if (length(sample_ids) >= 2) sample_ids[2] else first_sample

default_compare_left <- if ("BC01" %in% sample_ids) "BC01" else first_sample
default_compare_right <- if ("BC02" %in% sample_ids) "BC02" else second_sample

default_gene <- get_default_gene(first_sample, svg_dir)
if (is.null(default_gene)) default_gene <- NA_character_

ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    base_font = font_google("Space Grotesk"),
    heading_font = font_google("Space Grotesk"),
    code_font = font_google("JetBrains Mono"),
    primary = "#1f1f1f",
    bg = "#f6f6f6",
    fg = "#141414"
  ),
  tags$head(
    tags$title("Spatial Insights Cohort"),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  div(
    class = "app-shell",
    div(
      class = "app-hero",
      h2("Spatial Insights Cohort"),
      p("Business-ready overview of cohort spatial results")
    ),
    div(
      class = "panel-card control-card",
      div(class = "section-title", "Controls"),
      fluidRow(
        class = "control-row",
        column(4, selectInput("sample_id", "Sample", choices = sample_ids, selected = first_sample)),
        column(4, sliderInput("top_n", "Top N SVGs", min = 10, max = 50, value = 25, step = 1)),
        column(4, checkboxInput("advanced_svg", "Inspect single gene", value = FALSE))
      ),
      div(
        class = "control-actions",
        actionButton("reset_defaults", "Reset defaults")
      )
    ),
    tabsetPanel(
      tabPanel(
        "Overview",
        fluidRow(
          column(12, uiOutput("overview_kpis"))
        ),
        fluidRow(
          column(
            7,
            div(
              class = "panel-card",
              h4("Run Summary"),
              uiOutput("run_summary")
            )
          ),
          column(
            5,
            div(
              class = "panel-card",
              h4("Data readiness"),
              uiOutput("data_ready")
            )
          )
        ),
        div(
          class = "panel-card",
          h4("Summary"),
          p("This dashboard summarizes cohort outputs for quality checks, reporting, and stakeholder review.")
        ),
        bslib::accordion(
          bslib::accordion_panel(
            "Diagnostics",
            uiOutput("diagnostics_panel")
          ),
          open = NULL
        )
      ),
    tabPanel(
      "Domains",
      uiOutput("regions_summary"),
      p("Inferred spatial domains."),
      plotlyOutput("domains_plot", height = "600px"),
      bslib::accordion(
        bslib::accordion_panel(
          "Diagnostics",
          uiOutput("regions_troubleshoot")
        ),
        open = NULL
      )
    ),
    tabPanel(
      "Spatial Genes",
      uiOutput("svg_summary"),
      plotlyOutput("svg_score_plot", height = "600px"),
      bslib::accordion(
        bslib::accordion_panel(
          "Diagnostics",
          div(class = "panel-card", plotOutput("svg_top_bar", height = "200px"))
        ),
        open = NULL
      ),
      conditionalPanel(
        "input.advanced_svg == true",
        selectInput("gene", "Gene", choices = get_svg_genes(svg_dir, first_sample), selected = default_gene),
        plotlyOutput("svg_feature_plot", height = "500px")
      ),
      bslib::accordion(
        bslib::accordion_panel(
          "Troubleshooting",
          uiOutput("svg_troubleshoot")
        ),
        open = NULL
      )
    ),
    tabPanel(
      "Neighborhoods",
      uiOutput("nb_summary"),
      plotlyOutput("nb_heatmap", height = "600px"),
      bslib::accordion(
        bslib::accordion_panel(
          "Top domain pairs",
          DTOutput("nb_top_pairs")
        ),
        open = NULL
      ),
      bslib::accordion(
        bslib::accordion_panel(
          "Diagnostics",
          uiOutput("nb_troubleshoot")
        ),
        open = NULL
      )
    ),
    tabPanel(
      "Compare",
      div(
        class = "panel-card control-card",
        div(class = "section-title", "Domains comparison"),
        fluidRow(
          class = "control-row",
          column(6, selectInput("compare_left", "Left sample", choices = sample_ids, selected = default_compare_left)),
          column(6, selectInput("compare_right", "Right sample", choices = sample_ids, selected = default_compare_right))
        )
      ),
      fluidRow(
        column(6, plotlyOutput("compare_domains_left", height = "520px")),
        column(6, plotlyOutput("compare_domains_right", height = "520px"))
      )
    ),
    tabPanel(
      "Methods",
      tags$ul(
        tags$li("Data source: DLPFC Visium"),
        tags$li("QC: threshold-based spot filtering"),
        tags$li("Normalization: log-normalization + PCA"),
        tags$li("Domains: BayesSpace (if available) or k-means fallback"),
        tags$li("SVG: Moran's I on top variable genes"),
        tags$li("Neighborhood: kNN adjacency summaries")
      )
    )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    domain_frames = list(),
    svg_score_frames = list(),
    gene_frames = list(),
    nb_matrices = list(),
    nb_pairs = list()
  )

  observeEvent(input$sample_id, {
    sample_id <- input$sample_id
    def_gene <- get_default_gene(sample_id, svg_dir)
    genes <- get_svg_genes(svg_dir, sample_id)
    if (length(genes) > 0) {
      updateSelectInput(session, "gene", choices = genes, selected = def_gene)
    }
  }, ignoreInit = FALSE)

  observeEvent(input$reset_defaults, {
    updateSelectInput(session, "sample_id", selected = first_sample)
    updateSliderInput(session, "top_n", value = 25)
    updateCheckboxInput(session, "advanced_svg", value = FALSE)
    def_gene <- get_default_gene(first_sample, svg_dir)
    if (!is.null(def_gene)) updateSelectInput(session, "gene", selected = def_gene)
  })

  output$overview_kpis <- renderUI({
    cm <- compute_cohort_metrics()
    if (!isTRUE(cm$ok)) {
      return(div(class = "panel-card", "Cohort metrics not available."))
    }
    dom_median <- cm$metrics$domains_median
    div(
      class = "kpi-grid",
      div(class = "kpi-card", div(class = "kpi-label", "Samples"), div(class = "kpi-value", cm$metrics$n_samples)),
      div(class = "kpi-card", div(class = "kpi-label", "Domains per sample"), div(class = "kpi-value", dom_median)),
      div(class = "kpi-card", div(class = "kpi-label", "Mean stability"), div(class = "kpi-value", cm$metrics$mean_stability)),
      div(class = "kpi-card", div(class = "kpi-label", "Top domain pair"), div(class = "kpi-value", cm$metrics$top_consensus_pair))
    )
  })

  output$run_summary <- renderUI({
    cm <- compute_cohort_metrics()
    meta <- read_run_metadata(project_root)
    if (is.null(meta)) {
      return(div(
        tags$p(tags$strong("Run metadata not found.")),
        tags$p("Using cohort summaries as fallback."),
        if (isTRUE(cm$ok)) tags$ul(
          tags$li(paste("Mean stability:", cm$metrics$mean_stability)),
          tags$li(paste("Top recurring domain pair:", cm$metrics$top_consensus_pair))
        )
      ))
    }

    sample_ids <- meta$sample_ids
    samples_text <- "none"
    if (!is.null(sample_ids) && length(sample_ids) > 0) {
      if (length(sample_ids) <= 12) {
        samples_text <- paste(sample_ids, collapse = ", ")
      } else {
        samples_text <- paste(c(sample_ids[1:12], "and more"), collapse = ", ")
      }
    }

    run_date <- meta$run_timestamp
    steps_used <- meta$params$steps
    dataset_used <- meta$params$dataset

    div(
      tags$ul(
        tags$li(paste("Run date:", run_date)),
        tags$li(paste("Samples used:", samples_text)),
        tags$li(paste("Steps:", steps_used)),
        tags$li(paste("Dataset:", dataset_used)),
        tags$li(paste("Mean stability:", if (isTRUE(cm$ok)) cm$metrics$mean_stability else "NA")),
        tags$li(paste("Top recurring domain pair:", if (isTRUE(cm$ok)) cm$metrics$top_consensus_pair else "NA"))
      )
    )
  })

  output$data_ready <- renderUI({
    items <- list(
      Domains = file.exists(file.path(domains_dir, "domains_cohort_summary.csv")),
      "Spatial Genes" = file.exists(file.path(svg_dir, "svg_cohort_summary.csv")),
      Neighborhoods = file.exists(file.path(nb_dir, "neighborhood_cohort_summary.csv"))
    )
    tags$div(lapply(names(items), function(name) {
      ok <- items[[name]]
      tags$div(
        class = "status-row",
        tags$span(class = paste("status-dot", if (ok) "green" else "red")),
        tags$span(name)
      )
    }))
  })

  output$diagnostics_panel <- renderUI({
    sample_id <- input$sample_id
    qc_plot <- file.path(project_root, "outputs", "qc", "plots", paste0(sample_id, "_qc.png"))
    pca_plot <- file.path(project_root, "outputs", "integration", "plots", paste0(sample_id, "_pca_var.png"))
    items <- list(
      "QC histograms" = qc_plot,
      "PCA variance" = pca_plot
    )
    tags$ul(lapply(names(items), function(label) {
      path <- items[[label]]
      exists <- file.exists(path)
      tags$li(paste(label, if (exists) path else "not found"))
    }))
  })

  output$regions_summary <- renderUI({
    dm <- compute_domains_metrics(input$sample_id)
    if (!isTRUE(dm$ok)) {
      return(div(
        class = "panel-card",
        tags$strong("Missing domains data"),
        tags$p(dm$error)
      ))
    }
    div(
      class = "panel-card",
      h4(paste0("Domains: ", input$sample_id)),
      tags$ul(
        tags$li(paste("Domains:", dm$metrics$n_domains)),
        tags$li(paste("Largest region:", dm$metrics$largest_domain_pct, "%"))
      )
    )
  })

  output$domains_plot <- renderPlotly({
    sample_id <- input$sample_id
    withProgress(message = "Loading domains", value = 0.2, {
      if (!sample_id %in% names(rv$domain_frames)) {
        dom_path <- resolve_domains_rds(processed_dir, sample_id)
        if (!file.exists(dom_path)) {
          return(plot_placeholder_plotly("Domains not available"))
        }
        rv$domain_frames[[sample_id]] <- load_domain_frame(dom_path)
      }
      df <- rv$domain_frames[[sample_id]]
      if (is.null(df)) return(plot_placeholder_plotly("Domains not available"))
      gg <- plot_domain_scatter(df)
      plotly::ggplotly(gg, tooltip = "text")
    })
  })

  output$compare_domains_left <- renderPlotly({
    sample_id <- input$compare_left
    req(sample_id)
    withProgress(message = "Loading domains (left)", value = 0.2, {
      if (!sample_id %in% names(rv$domain_frames)) {
        dom_path <- resolve_domains_rds(processed_dir, sample_id)
        if (!file.exists(dom_path)) {
          return(plot_placeholder_plotly("Domains not available"))
        }
        rv$domain_frames[[sample_id]] <- load_domain_frame(dom_path)
      }
      df <- rv$domain_frames[[sample_id]]
      if (is.null(df)) return(plot_placeholder_plotly("Domains not available"))
      gg <- plot_domain_scatter(df)
      plotly::ggplotly(gg, tooltip = "text")
    })
  })

  output$compare_domains_right <- renderPlotly({
    sample_id <- input$compare_right
    req(sample_id)
    withProgress(message = "Loading domains (right)", value = 0.2, {
      if (!sample_id %in% names(rv$domain_frames)) {
        dom_path <- resolve_domains_rds(processed_dir, sample_id)
        if (!file.exists(dom_path)) {
          return(plot_placeholder_plotly("Domains not available"))
        }
        rv$domain_frames[[sample_id]] <- load_domain_frame(dom_path)
      }
      df <- rv$domain_frames[[sample_id]]
      if (is.null(df)) return(plot_placeholder_plotly("Domains not available"))
      gg <- plot_domain_scatter(df)
      plotly::ggplotly(gg, tooltip = "text")
    })
  })

  output$regions_troubleshoot <- renderUI({
    dom_path <- resolve_domains_rds(processed_dir, input$sample_id)
    div(
      tags$p(paste("Expected file:", dom_path))
    )
  })

  output$svg_summary <- renderUI({
    sm <- compute_svg_metrics(input$sample_id)
    if (!isTRUE(sm$ok)) {
      return(div(
        class = "panel-card",
        tags$strong("Missing SVG data"),
        tags$p(sm$error)
      ))
    }
    div(
      class = "panel-card",
      h4(paste0("Spatial genes: ", input$sample_id)),
      tags$ul(
        tags$li(paste("Top gene:", sm$metrics$top_gene)),
        tags$li(paste("Genes tested:", sm$metrics$n_genes_tested)),
        tags$li(paste("Consensus hit:", sm$metrics$consensus_hit))
      )
    )
  })

  output$svg_score_plot <- renderPlotly({
    sample_id <- input$sample_id
    top_n <- input$top_n
    withProgress(message = "Scoring SVG set", value = 0.2, {
      key <- paste0(sample_id, "_", top_n)
      if (!key %in% names(rv$svg_score_frames)) {
        df <- load_svg_score_frame(sample_id, processed_dir, svg_dir, top_n)
        if (is.null(df)) return(plot_placeholder_plotly("SVG score not available"))
        rv$svg_score_frames[[key]] <- df
      }
      df <- rv$svg_score_frames[[key]]
      gg <- plot_svg_score_scatter(df)
      plotly::ggplotly(gg, tooltip = "text")
    })
  })

  output$svg_top_bar <- renderPlot({
    df <- load_svg_ranked(input$sample_id, svg_dir)
    if (is.null(df)) return(invisible(NULL))
    plot_svg_bar(df, n = 10)
  })

  output$svg_feature_plot <- renderPlotly({
    if (!isTRUE(input$advanced_svg)) return(NULL)
    sample_id <- input$sample_id
    gene <- input$gene
    withProgress(message = "Loading gene", value = 0.2, {
      if (!sample_id %in% names(rv$gene_frames) || is.null(rv$gene_frames[[sample_id]][[gene]])) {
        df <- load_gene_values(sample_id, processed_dir, gene)
        if (is.null(df)) return(plot_placeholder_plotly("Gene data not available"))
        rv$gene_frames[[sample_id]][[gene]] <- df
      }
      df <- rv$gene_frames[[sample_id]][[gene]]
      gg <- plot_gene_scatter(df, gene)
      plotly::ggplotly(gg, tooltip = "text")
    })
  })

  output$svg_troubleshoot <- renderUI({
    svg_path <- file.path(svg_dir, paste0("svg_ranked_", input$sample_id, ".csv"))
    div(
      tags$p(paste("Expected file:", svg_path))
    )
  })

  output$nb_summary <- renderUI({
    nm <- compute_neighborhood_metrics(input$sample_id)
    if (!isTRUE(nm$ok)) {
      return(div(
        class = "panel-card",
        tags$strong("Missing neighborhood data"),
        tags$p(nm$error)
      ))
    }
    div(
      class = "panel-card",
      h4(paste0("Neighborhoods: ", input$sample_id)),
      tags$ul(
        tags$li(paste("Top pair:", nm$metrics$top_pair)),
        tags$li(paste("Top score:", nm$metrics$top_pair_score)),
        tags$li(paste("Matrix density:", nm$metrics$matrix_density))
      )
    )
  })

  output$nb_heatmap <- renderPlotly({
    sample_id <- input$sample_id
    withProgress(message = "Loading neighborhood", value = 0.2, {
      if (!sample_id %in% names(rv$nb_matrices)) {
        rv$nb_matrices[[sample_id]] <- load_neighborhood_matrix(sample_id, nb_dir)
      }
      mat <- rv$nb_matrices[[sample_id]]
      if (is.null(mat)) return(plot_placeholder_plotly("Neighborhood not available"))
      plot_heatmap_plotly(mat, paste0("Neighborhood adjacency: ", sample_id))
    })
  })

  output$nb_top_pairs <- renderDT({
    sample_id <- input$sample_id
    if (!sample_id %in% names(rv$nb_pairs)) {
      path <- file.path(nb_dir, paste0("neighborhood_top_pairs_", sample_id, ".csv"))
      rv$nb_pairs[[sample_id]] <- load_csv_safe(path)
    }
    df <- rv$nb_pairs[[sample_id]]
    if (is.null(df)) return(NULL)
    DT::datatable(head(df, 10), options = list(pageLength = 10))
  })

  output$nb_troubleshoot <- renderUI({
    nb_path <- file.path(nb_dir, paste0("neighborhood_matrix_", input$sample_id, ".csv"))
    div(
      tags$p(paste("Expected file:", nb_path))
    )
  })
}

shinyApp(ui, server)
