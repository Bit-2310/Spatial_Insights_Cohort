# UI Component Helpers for Shiny Dashboard

#' Create a KPI card
kpi_card <- function(label, value) {
    div(
        class = "kpi-card",
        div(class = "kpi-label", label),
        div(class = "kpi-value", value)
    )
}

#' Create a panel card wrapper
panel_card <- function(..., title = NULL) {
    content <- list(...)
    if (!is.null(title)) {
        content <- c(list(h4(title)), content)
    }
    div(class = "panel-card", content)
}

#' Create a status indicator
status_indicator <- function(label, available = TRUE) {
    div(
        class = "status-row",
        span(class = paste0("status-dot ", if (available) "green" else "red")),
        span(label)
    )
}

#' Create a placeholder plot for missing data
plot_placeholder <- function(message = "Data not available") {
    div(
        class = "plot-placeholder",
        style = "padding: 60px; text-align: center; color: #999;",
        tags$p(style = "font-size: 1.1em;", message)
    )
}
