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

#' Observation level simple model test.
#'
#' \code{alt1_smpl2} uses the results of a simple linear model in the probit
#' space estimated using the input set to estimate the observations described
#' by the test set. The differences are summarised by the mean squared error,
#' mean absolute deviation, and median absolute deviation.
#'
#' @param input Input set data frame.
#' @param tests Test set data frame.
#'
#' @family testing functions
#'
#' @return Returns a summary table with the model performance.
#'
#' @export
alt1_smpl2 <- function(input, tests) {
  if (nrow(tests) == 0) {
    empty <- input %>%
      select(variable, sex) %>%
      distinct %>%
      mutate(model = "ABC",
             mse.100 = NA,
             mae.100 = NA,
             mbe.100 = NA,
             mmad.100 = NA) %>%
      select(model, everything()) %>%
      rename(level = variable)
    return(empty)
  }

  df <- bind_rows(input %>% mutate(class = "train"),
                  tests %>% mutate(class = "test"))

  smpl <- function(x) {
    mod <- lm(data = x %>% filter(class == "train") %>% mutate(value = qnorm(value)), value ~ year)
    x %>%
      filter(class == "test") %>%
      mutate(value = qnorm(value)) %>%
      {mutate(., pred = predict(mod, type = 'response', newdata = .))} %>%
      mutate(value = pnorm(value), pred = pnorm(pred)) %>%
      select(-class)
  }

  df %>%
    group_by(country, variable, sex) %>%
    nest() %>%
    mutate(data = purrr::map(data, smpl)) %>%
    unnest(data) %>%
    ungroup() %>%
    mutate(diff = pred - value, diff_sq = diff^2) %>%
    group_by(variable, sex) %>%
    summarise(mse.100 = round(100 * mean(diff_sq, na.rm = TRUE), 3),
              mae.100 = round(100 * mean(abs(diff), na.rm = TRUE), 3),
              mbe.100 = round(100 * mean(diff, na.rm = TRUE), 3),
              mmad.100 = round(100 * median(abs(diff), na.rm = TRUE), 3)) %>%
    mutate(model = "simple", variable = variable, sex = sex) %>%
    select(model, level = variable, sex, everything())
}

#' Survey level simple model test.
#'
#' \code{alt1_smpl3} uses the results of a simple linear model in the probit
#' space estimated using the input set to estimate the completion rate indicator
#' of the observation year of each survey in the test set. The differences are
#' summarised by the mean squared error, mean absolute deviation, and median
#' absolute deviation.
#'
#' @param input Input set data frame.
#' @param tests Test set data frame.
#'
#' @family testing functions
#'
#' @return Returns a summary table with the model performance.
#'
#' @export
alt1_smpl3 <- function(input, tests) {
  if (nrow(tests) == 0) {
    empty <- input %>%
      select(variable, sex) %>%
      distinct %>%
      mutate(model = "ABC",
             mse.100 = NA,
             mae.100 = NA,
             mbe.100 = NA,
             mmad.100 = NA) %>%
      select(model, everything()) %>%
      rename(level = variable)
    return(empty)
  }

  df <- bind_rows(input %>% mutate(class = "train"),
                  tests %>% mutate(class = "test"))

  smpl <- function(x) {
    mod <- lm(data = x %>% filter(class == "train") %>% mutate(value = qnorm(value)), value ~ year)
    x %>%
      filter(class == "test") %>%
      mutate(value = qnorm(value)) %>%
      {mutate(., pred = predict(mod, type = 'response', newdata = .))} %>%
      mutate(value = pnorm(value), pred = pnorm(pred)) %>%
      select(-class)
  }

  df %>%
    group_by(country, variable, sex) %>%
    nest() %>%
    mutate(data = purrr::map(data, smpl)) %>%
    unnest(data) %>%
    ungroup() %>%
    group_by(country, variable, sex) %>%
    filter(recondist == min(recondist)) %>%
    summarise(pred = mean(pred), value = mean(value)) %>%
    mutate(diff = pred - value, diff_sq = diff^2) %>%
    group_by(variable, sex) %>%
    summarise(mse.100 = round(100 * mean(diff_sq, na.rm = TRUE), 3),
              mae.100 = round(100 * mean(abs(diff), na.rm = TRUE), 3),
              mbe.100 = round(100 * mean(diff, na.rm = TRUE), 3),
              mmad.100 = round(100 * median(abs(diff), na.rm = TRUE), 3)) %>%
    mutate(model = "simple", variable = variable, sex = sex) %>%
    select(model, level = variable, sex, everything())
}

#' Observation level flat model test.
#'
#' \code{alt2_flat2} uses the results of an intercept only model in the probit
#' space estimated using the input set to estimate the observations described
#' by the test set. The differences are summarised by the mean squared error,
#' mean absolute deviation, and median absolute deviation.
#'
#' @param input Input set data frame.
#' @param tests Test set data frame.
#'
#' @family testing functions
#'
#' @return Returns a summary table with the model performance.
#'
#' @export
alt2_flat2 <- function(input, tests) {
  if (nrow(tests) == 0) {
    empty <- input %>%
      select(variable, sex) %>%
      distinct %>%
      mutate(model = "ABC",
             mse.100 = NA,
             mae.100 = NA,
             mbe.100 = NA,
             mmad.100 = NA) %>%
      select(model, everything()) %>%
      rename(level = variable)
    return(empty)
  }

  df <- bind_rows(input %>% mutate(class = "train"),
                  tests %>% mutate(class = "test"))

  flat <- function(x) {
    mod <- lm(data = x %>% filter(class == "train") %>% mutate(value = qnorm(value)), value ~ 1)
    x %>%
      filter(class == "test") %>%
      mutate(value = qnorm(value)) %>%
      {mutate(., pred = predict(mod, type = 'response', newdata = .))} %>%
      mutate(value = pnorm(value), pred = pnorm(pred)) %>%
      select(-class)
  }

  df %>%
    group_by(country, variable, sex) %>%
    nest() %>%
    mutate(data = purrr::map(data, flat)) %>%
    unnest(data) %>%
    ungroup() %>%
    mutate(diff = pred - value, diff_sq = diff^2) %>%
    group_by(variable, sex) %>%
    summarise(mse.100 = round(100 * mean(diff_sq, na.rm = TRUE), 3),
              mae.100 = round(100 * mean(abs(diff), na.rm = TRUE), 3),
              mbe.100 = round(100 * mean(diff, na.rm = TRUE), 3),
              mmad.100 = round(100 * median(abs(diff), na.rm = TRUE), 3)) %>%
    mutate(model = "flat", variable = variable, sex = sex) %>%
    select(model, level = variable, sex, everything())
}

#' Survey level flat model test.
#'
#' \code{alt2_flat3} uses the results of an intercept only model in the probit
#' space estimated using the input set to estimate the completion rate indicator
#' of the observation year of each survey in the test set. The differences are
#' summarised by the mean squared error, mean absolute deviation, and median
#' absolute deviation.
#'
#' @param input Input set data frame.
#' @param tests Test set data frame.
#'
#' @family testing functions
#'
#' @return Returns a summary table with the model performance.
#'
#' @export
alt2_flat3 <- function(input, tests) {
  if (nrow(tests) == 0) {
    empty <- input %>%
      select(variable, sex) %>%
      distinct %>%
      mutate(model = "ABC",
             mse.100 = NA,
             mae.100 = NA,
             mbe.100 = NA,
             mmad.100 = NA) %>%
      select(model, everything()) %>%
      rename(level = variable)
    return(empty)
  }

  df <- bind_rows(input %>% mutate(class = "train"),
                  tests %>% mutate(class = "test"))

  flat <- function(x) {
    mod <- lm(data = x %>% filter(class == "train") %>% mutate(value = qnorm(value)), value ~ 1)
    x %>%
      filter(class == "test") %>%
      mutate(value = qnorm(value)) %>%
      {mutate(., pred = predict(mod, type = 'response', newdata = .))} %>%
      mutate(value = pnorm(value), pred = pnorm(pred)) %>%
      select(-class)
  }

  df %>%
    group_by(country, variable, sex) %>%
    nest() %>%
    mutate(data = purrr::map(data, flat)) %>%
    unnest(data) %>%
    ungroup() %>%
    group_by(country, variable, sex) %>%
    filter(recondist == min(recondist)) %>%
    summarise(pred = mean(pred), value = mean(value)) %>%
    mutate(diff = pred - value, diff_sq = diff^2) %>%
    group_by(variable, sex) %>%
    summarise(mse.100 = round(100 * mean(diff_sq, na.rm = TRUE), 3),
              mae.100 = round(100 * mean(abs(diff), na.rm = TRUE), 3),
              mbe.100 = round(100 * mean(diff, na.rm = TRUE), 3),
              mmad.100 = round(100 * median(abs(diff), na.rm = TRUE), 3)) %>%
    mutate(model = "flat", variable = variable, sex = sex) %>%
    select(model, level = variable, sex, everything())
}

#' Survey level latest model test.
#'
#' \code{alt3_ltst3} uses the results of a "latest available value" model
#' estimated using the input set to estimate the completion rate indicator
#' of the observation year of each survey in the test set. The differences
#' are summarised by the mean squared error,mean absolute deviation, and
#' median absolute deviation.
#'
#' @param input Input set data frame.
#' @param tests Test set data frame.
#'
#' @family testing functions
#'
#' @return Returns a summary table with the model performance.
#'
#' @export
alt3_ltst3 <- function(input, tests) {
  preds <- input %>%
    group_by(country, variable, sex) %>%
    filter(!is.na(value)) %>%
    filter(obsyear == max(obsyear) & recondist == min(recondist)) %>%
    summarise(pred = mean(value, na.rm = TRUE)) %>%
    ungroup

  tests %>%
    group_by(country, variable, sex) %>%
    filter(!is.na(value)) %>%
    filter(recondist == min(recondist)) %>%
    summarise(value = mean(value)) %>%
    left_join(preds, by = c("country", "variable", "sex")) %>%
    mutate(diff = pred - value, diff_sq = diff^2) %>%
    group_by(variable, sex) %>%
    summarise(mse.100 = round(100 * mean(diff_sq, na.rm = TRUE), 3),
              mae.100 = round(100 * mean(abs(diff), na.rm = TRUE), 3),
              mbe.100 = round(100 * mean(diff, na.rm = TRUE), 3),
              mmad.100 = round(100 * median(abs(diff), na.rm = TRUE), 3)) %>%
    mutate(model = "latest", variable = variable, sex = sex) %>%
    select(model, level = variable, sex, everything())
}
