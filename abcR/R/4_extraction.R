#' Extract observation replications.
#'
#' @param input Input data frame.
#' @param mcmc ABC model output.
#'
#' @return Data frame containing the replications generated at each iteration
#' of the ABC model.
#'
#' @family extraction functions
#'
#' @export
get_yrepl <- function(input, mcmc) {
  input %>%
    mutate(index = row_number()) %>%
    select(index, country, year, survey) %>%
    inner_join(mcmc %>%
                 tidybayes::gather_draws(yrepl[n]) %>%
                 ungroup %>%
                 select(index=n, iteration=.draw, value=.value),
               by = c("index"))
}

#' Extract long-term drift
#'
#' Extracts drift iterations from the model output and summarises the iterations
#' in a 90\% interval.
#'
#' @inheritParams get_muzero
#'
#' @return Data frame with columns for country, value, and lower/upper bounds. The
#'   output does not specify the level and sex combination as it is assumed to be
#'   known by the user based on the inputs.
#'
#' @family extraction functions
#'
#' @export
get_drift <- function(df, raw_mcmc) {
  raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::gather_draws(drift[country]) %>%
    tidybayes::point_interval(.value, .width = 0.9) %>%
    #mutate(variable = level, sex = sex) %>%
    select(country, .value, .lower, .upper) %>%
    rename(value = .value, lower = .lower, upper = .upper) %>%
    ungroup
}

#' Extract Latent Space Intercept
#'
#' Extracts the intercept term from the model output and summarises the iterations
#' in a 90\% interval. Note that the intercept for the ABC dataset corresponds to
#' 1980.
#'
#' @param df Input data frame.
#' @param raw_mcmc ABC model output.
#'
#' @return Data frame with columns for country, value, and lower/upper bounds. The
#'   output does not specify the level and sex combination as it is assumed to be
#'   known by the user based on the inputs.
#'
#' @family extraction functions
#'
#' @export
get_muzero <- function(df, raw_mcmc) {
  raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::gather_draws(muzero[country]) %>%
    tidybayes::point_interval(.value, .width = 0.9) %>%
    select(country, .value, .lower, .upper) %>%
    rename(value = .value, lower = .lower, upper = .upper) %>%
    ungroup
}

#' Extract Survey Bias
#'
#' Extracts survey bias iterations from the model output and attaches information
#' about the survey. Note that no summary statistics are computed.
#'
#' @inheritParams get_muzero
#' @param rounds Survey rounds data frame.
#'
#' @return Data frame with columns for survey, iteration, variable, value, country,
#'   year, round, level, and sex.
#'
#' @family extraction functions
#'
#' @export
get_bias <- function(df, raw_mcmc, rounds) {
  raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::gather_draws(beta_s[survey]) %>%
    ungroup %>%
    select(-.iteration, -.chain) %>%
    rename(iteration = .draw, variable = .variable, value = .value) %>%
    left_join(mutate(rounds, survey = paste0(country, survey, year))) %>%
    mutate(level = (df %>% ungroup %>% select(variable) %>% distinct %>% pull),
           sex = (df %>% ungroup %>% select(sex) %>% distinct %>% pull))
}

#' Extract Model Parameters
#'
#' Extracts survey, late, very late, and age misreporting bias terms from the
#' model output, then presents the median and 90 percent intervals for each parameter.
#'
#' @param df Input data frame.
#' @param raw_mcmc ABC model output.
#'
#' @return Data frame with columns for survey, parameter, country, value, lower, upper,
#'   and level.
#'
#' @family extraction functions
#'
#' @export
get_pars <- function(df, raw_mcmc) {
  raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::gather_draws(
      beta_s[survey],
      late[country],
      mult5err[country],
      vlate[country]
    ) %>%
    ungroup %>%
    select(-.iteration, -.chain) %>%
    rename(iteration = .draw, variable = .variable, value = .value) %>%
    group_by(survey, variable, country) %>%
    tidybayes::point_interval(value, .width = 0.9) %>%
    select(-.width, -.point, -.interval) %>%
    rename(parameter = variable, lower = .lower, upper = .upper) %>%
    ungroup %>%
    mutate(level = (df %>% ungroup %>% select(variable) %>% distinct %>% pull))
}

#' Extract Rhats
#'
#' Extracts Rhats for all parameters (base and transformed). Note that generated quantities
#' are excluded. The 9 parameters with the worst Rhats are returned.
#'
#' @param df Input data frame.
#' @param level Corresponding level.
#' @param sex Corresponding sex.
#'
#' @family extraction functions
#'
#' @export
get_rhats <- function(df, level, sex) {
  rhat      <- bayesplot::rhat(df)
  top9_rhat <- utils::tail(sort(rhat[!stringr::str_detect(names(rhat), c('yhat|yrepl|ll'))]), 9)

  #tibble(parameter = ordered(names(top9_rhat)), rhat = top9_rhat, level = level, sex = sex) %>%
  #  mutate(parameter = reorder(parameter, rhat))
  tibble(parameter = names(top9_rhat), rhat = top9_rhat, level = level, sex = sex)
}

#' Extract Parameter Samples.
#'
#' Extracts samples for a specified parameter family to flow into the generated quantities
#' model.
#'
#' @param raw_mcmc Stan model output
#' @param par Parameter name
#'
#' @family extraction functions
get_parsamps <- function(raw_mcmc, par) {
  extract(mcmc)[["par"]]
}

#' Extract Non-Sampling Variance Samples.
#'
#' Extracts samples for non-sampling variance and subsequently merges them by country
#' for input into the leave-survey-out generated quantities model.
#'
#' @param df Input data frame
#' @param raw_mcmc Stan model output
#'
#' @family extraction functions
get_nsvar <- function(df, raw_mcmc) {
  unname(as.matrix(raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::gather_draws(
      sigma_s[survey]
    ) %>%
    ungroup %>%
    mutate(country = substr(survey, 1, 3)) %>%
    select(-.chain, -.iteration, -.variable, -survey) %>%
    rename(iteration = .draw, value = .value) %>%
    group_by(country, iteration) %>%
    summarise(m = mean(value)) %>%
    pivot_wider(names_from = "country", values_from = "m") %>%
    select(-iteration) %>%
    select(order(colnames(.)))))
}
