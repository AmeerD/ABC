#' abcR: A package for running the Adjusted Bayesian Completion Rates
#' (ABC) model pipeline.
#'
#' The abcR package provides categories of functions for each stage
#' of the Adjusted Bayesian Completion Rates (ABC) model:
#' \enumerate{
#'   \item Data compilation
#'   \item Data pre-processing
#'   \item Modeling
#'   \item Post-processing
#'   \item Output Extraction
#'   \item Visualization
#' }
#'
#' As a general rule, examples are not shown in the documentation for
#' this package. The reason being that the vast majority of the
#' functions are written to operate on specific dataframes or Stan results
#' and their outputs chain together. These dataframes are in most cases far
#' too large to share within the package. Instead, we refer the user to the
#' following set of \href{https://github.com/ropensci/drake}{drake} plans that
#' illustrate how to use the package: INSERT LINK TO PLAN FOLDER HERE
#'
#' @section Output extraction:
#' The function for output extraction are:
#' \enumerate{
#'   \item \code{\link{get_muzero}}
#'   \item \code{\link{get_bias}}
#' }
#'
#' @docType package
#' @import dplyr
#' @import ggplot2
#' @importFrom stats median dnorm qnorm pnorm rnorm runif rexp
#' @importFrom stats lm predict sd quantile weighted.mean na.omit
#' @importFrom tidyr nest unnest
#' @importFrom rstan stan_model stan sampling
#' @name abcR
NULL
