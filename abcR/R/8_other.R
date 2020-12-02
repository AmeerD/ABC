#' Plot walker?
#'
#' Saves plots for each country. Replaced by report drake target?
#'
#' @param df Output data frame.
#' @param pl_fun Plotting function.
#' @param lbl File name label.
#' @param target_yr Final year to plot.
#'
#' @export
plot_walker <- function(df, pl_fun, lbl, target_yr = 2020) {
  purrr::walk(unique(df$observed$country),
              ~ {print(.x); ggsave(plot = pl_fun(df, .x, target_yr = target_yr),
                                   filename = paste0(lbl, '/', .x, '.png'), dpi = 100); print(.x)})
}

#' GPI plot?
#'
#' @param d Output data frame.
#' @param c Countries to plot.
#'
#' @export
gpi_plot <- function(d, c) {
  mdl_plot(d, c)+
    scale_y_continuous(limits = c(0, 2))+
    geom_hline(yintercept = 1, colour = 'red')
}
