#' Censor latest (Deprecated).
#'
#' This function sets the value of observations attributed to the most recent
#' survey for each country to NA. These censored values are used as a test set
#' with which the ABC model can be evaluated.
#'
#' @param df Input data frame.
#'
#' @return Returns df with each test set observation having NA as its value.
#' @export
cnsr_ltst <- function(df) {
  df %>%
    group_by(country) %>%
    # filter(n_distinct(obsyear) > 2) %>%
    mutate(value = ifelse(obsyear == max(obsyear), NA, value)) %>%
    ungroup
}

#' Simple model test (Deprecated).
#'
#' Executes a simple linear model in the probit space.
#'
#' @param df Input data frame.
#'
#' @export
alt1_smpl <- function(df) {
  smpl <- function(x) {
    x %>%
      mutate(value = qnorm(value)) %>%
      {mutate(., pred = predict(lm(data = ., value ~ year), type = 'response', newdata = .))} %>%
      mutate(value = pnorm(value), pred = pnorm(pred))
  }

  df %>%
    nest(-country, -variable) %>%
    mutate(data = purrr::map(data, smpl)) %>%
    unnest %>%
    mutate(model = 'simple') %>%
    select(country, level = variable, model, year, value, pred) %>%
    arrange(year)
}

#' Flat model test (Deprecated).
#'
#' Executes an intercept only model in the probit space.
#'
#' @param df Input data frame.
#'
#' @export
alt2_flat <- function(df) {
  flat <- function(x) {
    x %>%
      mutate(value = qnorm(value)) %>%
      {mutate(., pred = predict(lm(data = ., value ~ 1), type = 'response', newdata = .))} %>%
      mutate(value = pnorm(value), pred = pnorm(pred))
  }

  df %>%
    nest(-country, -variable) %>%
    mutate(data = purrr::map(data, flat)) %>%
    unnest %>%
    mutate(model = 'flat') %>%
    select(country, level = variable, model, year, value, pred) %>%
    arrange(year)
}

#' Latest model test (Deprecated).
#'
#' Executes a latest survey model in the probit space.
#'
#' @param df Input data frame.
#'
#' @export
alt3_ltst <- function(df) {
  df %>%
    group_by(country, variable) %>%
    filter(!is.na(value)) %>%
    filter(
      obsyear == max(obsyear) &
        -0.5 < recondist & recondist < 2.5) %>%
    summarise(pred = mean(value, na.rm = TRUE)) %>%
    mutate(model = 'latest') %>%
    left_join(distinct(filter(df, is.na(value)), country, variable, survey)) %>%
    select(country, level = variable, model, pred)
}

# tp <- function(c) {
#   smpl_plot(err1, c)+geom_point(data = filter(df_errtgt, country == c), aes(y = value), shape = 15)
# }

#' Prediction summary (Deprecated).
#'
#' @param dfs List of prediction data frames.
#' @param errtgt Testing target.
#'
#' @export
pred_mse <- function(dfs, errtgt) {
  purrr::map_dfr(dfs, ~ right_join(.x, errtgt)) %>%
    group_by(model, level) %>%
    summarise(
      mse = round(1000 * mean((pred - value)^2, na.rm = TRUE), 1),
      mad = round(100 * mean(abs(pred - value), na.rm = TRUE), 1),
      mmad = round(100 * median(abs(pred - value), na.rm = TRUE), 1))
}

#' Observation level prediction summary (Deprecated).
#'
#' @param dfs List of prediction data frames.
#' @param errtgt Testing target.
#'
#' @export
pred_mse_ace <- function(dfs, errtgt) {
  target <-errtgt

  pred_mse(
    purrr::map(dfs, ~ select(., -contains('value'))),
    target) %>%
    rename(mse_ace = mse, mad_ace = mad, mmad_ace = mmad) %>%
    mutate(level = ordered(level, levels = c('usec', 'lsec', 'prim'),
                           labels = c('upper secondary', 'lower secondary', 'primary'))) %>%
    arrange(model, level)
}

#' Survey level prediction summary (Deprecated).
#'
#' @param dfs List of prediction data frames.
#' @param errtgt Testing target.
#'
#' @export
pred_mse_cre <- function(dfs, errtgt) {
  target <-
    errtgt %>%
    filter(-0.5 < obsage & obsage < 2.5)

  pred_mse(
    purrr::map(dfs, ~ {select(., -matches('value')) %>%
        semi_join(target) %>%
        group_by(model, level, country) %>%
        summarise(pred = mean(pred, na.rm = TRUE)) %>%
        ungroup})
    ,
    target %>%
      group_by(level, country) %>%
      summarise(value = mean(value, na.rm = TRUE)) %>%
      ungroup
  ) %>%
    rename(mse_cre = mse, mad_cre = mad) %>%
    mutate(level = ordered(level, levels = c('usec', 'lsec', 'prim'),
                           labels = c('upper secondary', 'lower secondary', 'primary'))) %>%
    arrange(model, level)
}
