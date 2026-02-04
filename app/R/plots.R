plot_domain_scatter <- function(df) {
  if (is.null(df)) return(NULL)
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = domain, text = paste0("barcode: ", barcode, "<br>domain: ", domain, "<br>sample: ", sample_id))) +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "x", y = "y", color = "domain")
}

plot_gene_scatter <- function(df, gene) {
  if (is.null(df)) return(NULL)
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = value, text = paste0("barcode: ", barcode, "<br>", gene, ": ", round(value, 3)))) +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "x", y = "y", color = gene)
}

plot_svg_score_scatter <- function(df) {
  if (is.null(df)) return(NULL)
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = score, text = paste0("barcode: ", barcode, "<br>score: ", round(score, 3)))) +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::geom_point(data = df[df$hotspot, ], color = "#000000", size = 0.8, shape = 21, stroke = 0.6) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "x", y = "y", color = "SVG score")
}

plot_svg_bar <- function(df, n = 10) {
  if (is.null(df)) return(NULL)
  score_col <- if ("score" %in% colnames(df)) "score" else if ("moran_I" %in% colnames(df)) "moran_I" else NULL
  if (is.null(score_col)) return(NULL)
  top <- head(df, n)
  ggplot2::ggplot(top, ggplot2::aes(x = reorder(gene, .data[[score_col]]), y = .data[[score_col]])) +
    ggplot2::geom_col(fill = "#0f6f6f") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "", y = "Score")
}

plot_heatmap_plotly <- function(mat, title) {
  if (is.null(mat)) return(NULL)
  plotly::plot_ly(
    x = colnames(mat),
    y = rownames(mat),
    z = mat,
    type = "heatmap",
    colorscale = "Viridis"
  ) |>
    plotly::layout(title = title, xaxis = list(title = ""), yaxis = list(title = ""))
}

plot_placeholder_plotly <- function(title, subtitle = NULL) {
  msg <- if (is.null(subtitle)) title else paste0(title, "<br>", subtitle)
  plotly::plot_ly() |>
    plotly::layout(
      title = list(text = msg),
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE)
    )
}
