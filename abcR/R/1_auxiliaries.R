
# AUXILIARIES -------------------------------------------------------------

#' Nonempty column check
#
#' @param c Data frame column.
#' @return Boolean.
#' @export
nonempty_col <- function(c) !all(is.na(c))

#' Add ordering to education levels
#'
#' @param df Input data frame.
#'
#' @export
order_levels <- function(df) {
  mutate(df, level = ordered(level, levels = c('upper secondary', 'lower secondary', 'primary')))
}

#' Substitute levels with full titles.
#'
#' @param ... Pass through parameters.
#'
#' @export
level_names <- function(...) {
  args <- list(...)
  names(args) <-
    all.vars(substitute(list(...))) %>%
    {case_when(
      str_detect(., 'usec') ~ 'upper secondary',
      str_detect(., 'lsec') ~ 'lower secondary',
      str_detect(., 'prim') ~ 'primary',
      TRUE ~ .)}
  args
}

# AGGREGATION AND VARIANCE ESTIMATION -------------------------------------

#' Aggregate surveys
#'
#' @param survey Survey data file.
#' @param bycats Category for aggregation.
#' @param estimate_se Boolean indicating whether to estimate standard errors.
#'
#' @export
agg_survey <- function(survey, bycats, estimate_se = FALSE) {
  path = paste0('data/', survey, '_prepped.qs');
  qs::qread(path) %>%
    filter(between(age, 10, 45)) %>%
    mutate(variable = level) %>%
    wide_aggregate(categories = bycats, depth = length(bycats), estimate_se = estimate_se) %>%
    filter(category == paste(bycats, collapse = ' and ')) %>%
    select(-category)
}

#' Load raw data into R.
#'
#' @param o_name Object name.
#' @param d_name Data name.
#'
#' @export
load_raw <- function(o_name, d_name) {suppressWarnings({
  rda_target <- glue::glue('data/{d_name}.rda')
  try(load(rda_target), silent = TRUE)
  if (!exists(o_name)) {
    try(silent = TRUE, {
      assign(o_name, haven::read_dta(glue::glue('data/{d_name}.dta')), envir = .GlobalEnv)
      save(o_name, file = rda_target)
    })}
  if (exists(o_name)) {
    print(glue::glue('Raw data successfully loaded into object {o_name}!'))
  } else {
    print(glue::glue('ERROR: No raw data {d_name}.rda or .dta available!'))
  }
})
}

#' Jackknife variance estimation
#'
#' @param d Input data frame.
#' @param estimate_se Boolean indicating whether to estimate standard errors.
#'
#' @export
wide_jk <- function(d, estimate_se = FALSE) {

  fun <- function(x) {(0.5 + sum(x$value * x$weight))/(1 + sum(x$weight)) #Beta(0.5,0.5) prior
    #weighted.mean(x$value, x$weight, na.rm = TRUE)
    }
  within(data.frame(value = fun(d)), {
    if (estimate_se) {
      # jk       <- jackknife(d)
      jk       <- modelr::crossv_loo(d)
      k        <- nrow(jk)
      # jk_ests  <- map_dbl(jk$sample, ~fun(as.data.frame(.)))
      jk_ests  <- purrr::map_dbl(jk$train, ~fun(as.data.frame(.)))
      se   <- sqrt(sum(((k-1) * value - (k - 1) * jk_ests)^2) / (k * (k - 1)))
      alph <- 0.5 + sum(d$value * d$weight)
      beta <- 1 + sum(d$weight)
      se_b <- sqrt((alph*beta)/((alph + beta)^2 * (alph + beta + 1)))
      #se_q <- sqrt(sum(((k-1) * qnorm(value) - (k - 1) * qnorm(jk_ests))^2) / (k * (k - 1)))
      rm(jk, k, jk_ests, alph, beta)
    }})
}

#' Aggregate for a given category set
#'
#' @param df Input data frame.
#' @param cs Categories.
#' @param estimate_se Boolean indicating whether to estimate standard errors.
#'
#' @export
wide_aggregate_by_cats <- function(df, cs, estimate_se = FALSE) {
  if (identical(cs, c(""))) {
    cats     <- NULL
    category <- "Total"
    n_dims   <- 0
  } else {
    cats     <- syms(cs)
    category <- paste(cats, collapse = ' and ')
    n_dims   <- length(cs)
  }

  print(category)

  df %>%
    group_by(!!! cats, country, year, survey, round, variable, cluster) %>%
    summarise(
      value = weighted.mean(value, weight, na.rm = TRUE),
      weight = sum(weight, na.rm = TRUE)
    ) %>%
    nest %>%
    mutate(data = purrr::map(data, ~ wide_jk(., estimate_se))) %>%
    unnest %>%
    na.omit %>%
    mutate(category = category, n_dims = n_dims) %>%
    select(country, year, survey, round, variable, value, category, everything()) %>%
    identity
}

#' Aggregate Over All Category Sets
#'
#' @param df Input data frame.
#' @param categories Aggregation categories.
#' @param depth Tree depth.
#' @param estimate_se Boolean indicating whether to estimate standard errors.
#'
#' @export
wide_aggregate <- function(df, categories = "", depth = 1, estimate_se = FALSE) {
  disaggs <- unique(c('', purrr::lmap(1:depth, function(n) utils::combn(categories, n, simplify = FALSE))))
  purrr::map_dfr(disaggs, function(c) wide_aggregate_by_cats(df, c, estimate_se))
}

#' Survey whitelist.
#'
#' Removes identified bad survey inputs.
#'
#' @param df Input data frame.
#'
#' @export
whitelist <- function(df) {
  filter(df,
          # PATCH: drop youngest from these presumed erroneous surveys
          !(survey %in% c('CANCensus2011',
                          'CHECensus2000',
                          'FRACensus2011',
                          'FRACensus2006',
                          'FRASILC2013',
                          'HUNCensus2001',
                          'ROUCensus2011',
                          'BLRCensus2009') & recondist <= 0),
          !(variable == 'prim' & (
            (survey %in% c('CANCensus2011', 'SDNCensus2008', 'SSDMICS2000', 'URYCensus2006',
                           'FRACensus2006', 'FRACensus2011', 'LTUSILC2013', 'BGRSILC2013',
                           'BELSILC2013', 'LVASILC2005'
                           # ,'MLTSILC2009' ##NEW ADDITION
            )) |
              (survey %in% c('LVASILC2005', 'MNGMICS2000', 'POLCensus2002', 'CHNCensus2000',
                             'FRASILC2017') & obsage > 0) |
              (survey %in% c('TTOMICS2011') & recondist > 10)
          )),
          !(variable == 'lsec' & ((survey %in% c('ALBDHS2009', 'ARMDHS2000', 'CANCensus2011', 'CHESILC2009',
                                                 'KENCensus2009', 'LVASILC2005', 'SDNCensus2008', 'SSDMICS2000', 'SVNSILC2005',
                                                 'VENMICS2000')) |
                                    (survey %in% c('HUNCensus2011', 'MNGMICS2000', 'POLCensus2002', 'ZWECensus2012') & obsage > 0) |
                                    (survey %in% c('BIHMICS2011') & recondist > 7)
          )),
          !(variable == 'usec' & ((survey %in% c('ARMDHS2016', 'ARMDHS2000', 'BIHMICS2011', 'EGYDHS2005', 'EGYDHS2014',
                                                 'KAZMICS2006', 'MDAMICS2000', 'NORSILC2009', 'NORSILC2013', 'POLCensus2002',
                                                 'POLCensus2011', 'SDNCensus2008', 'SRBSILC2013', 'SSDMICS2000', 'UKRCensus2001',
                                                 'KGZCensus2009', 'KAZMICS2006', 'BLRCensus2009',
                                                 'CHLCASEN2000', 'CHLCASEN2011', 'CHLCASEN2013', 'CHLCASEN2015')) |
                                    (survey %in% c('BELSILC2017') & obsage > 0) |
                                    (survey %in% c('KENDHS2003', 'TJKMICS2005') & recondist > 7) |
                                    (survey %in% c('TKMMICS2006', 'TKMMICS2015') & year >= 2002 & year <= 2014)
          ))
)}

#' Survey blacklist.
#'
#' @param df Input data frame
#'
#' @export
blacklist <- function(df) {
  anti_join(df, whitelist(df))
}

