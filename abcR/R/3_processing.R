#' Process ABC model output.
#'
#' \code{mdl_process} operates on the raw stanfit object outputted by the ABC
#' model and summarises the findings in a digestible format. Specifically,
#' variants of the latent true completion rate are computed and presented in
#' 90\% intervals.
#'
#' @param df Input data frame.
#' @param raw_mcmc Raw ABC model MCMC output.
#'
#' @return Data frame with columns for country, year, series, value, lower,
#' upper, variable, and sex. Note that series corresponds to one of:
#' \enumerate{
#'     \item projected3t5 (Official completion rate indicator)
#'     \item projected5 (Upper bracket of the indicator range)
#'     \item projected8 (Ultimate completion)
#' }
#'
#' @export
mdl_process <- function(df, raw_mcmc) {
  startyear <- min(df$year)
  baseyear <- startyear - 1
  variable <- df %>% ungroup %>% select(variable) %>% distinct() %>% pull
  sex <- df %>% ungroup %>% select(sex) %>% distinct() %>% pull

  lookup <- df %>%
    select(country, cap_adj) %>%
    distinct

  projected <- raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::spread_draws(
      mu_ct[country, year],
      late[country],
      vlate[country]
    ) %>%
    ungroup %>%
    mutate(year = year + baseyear) %>%
    select(-.iteration, -.chain) %>%
    rename(iteration = .draw) %>%
    group_by(country, iteration) %>%
    arrange(year) %>%
    mutate(
      mu5.ct = mu_ct,
      mu4.ct = dplyr::lead(mu_ct, 1) - late,
      mu3.ct = dplyr::lead(mu_ct, 2) - 2*late,
      mu8.ct = mu_ct + vlate * 3
    ) %>%
    select(-mu_ct, -late, -vlate) %>%
    mutate_at(.vars = vars(mu5.ct, mu4.ct, mu3.ct, mu8.ct), pnorm) %>%
    left_join(lookup, by="country") %>%
    mutate(mu5.ct = pmax(pmin(mu5.ct + cap_adj, 1), 0),
           mu4.ct = pmax(pmin(mu4.ct + cap_adj, 1), 0),
           mu3.ct = pmax(pmin(mu3.ct + cap_adj, 1), 0),
           mu8.ct = pmax(pmin(mu8.ct + cap_adj, 1), 0)) %>%
    select(-cap_adj) %>%
    mutate(projected3t5 = (mu3.ct + mu4.ct + mu5.ct)/3) %>%
    rename(projected5 = mu5.ct, projected8 = mu8.ct) %>%
    select(country, year, iteration, projected3t5, projected5, projected8) %>%
    tidyr::gather(series, value, projected3t5, projected5, projected8) %>%
    group_by(country, year, series) %>%
    tidybayes::point_interval(value, .width = 0.9) %>%
    select(-.width, -.point, -.interval) %>%
    rename(lower = .lower, upper = .upper) %>%
    ungroup %>%
    mutate_at(vars(value, lower, upper), ~ round(., 3)) %>%
    filter(year >= startyear)

  projected %>% mutate(variable = variable, sex = sex)
}

#' Process ABC model output for dropout model.
#'
#' \code{mdl_process_dropout} operates on the raw stanfit object outputted by the ABC
#' model and summarises the findings in a digestible format. Specifically,
#' variants of the latent true completion rate are computed and presented in
#' 90\% intervals.
#'
#' @param df Input data frame.
#' @param raw_mcmc Raw ABC model MCMC output.
#'
#' @return Data frame with columns for country, year, series, value, lower,
#' upper, variable, and sex. Note that series corresponds to one of:
#'
#' @export
mdl_process_dropout <- function(df, raw_mcmc) {
  startyear <- min(df$year)
  baseyear <- startyear - 1
  variable <- df %>% ungroup %>% select(variable) %>% distinct() %>% pull
  sex <- df %>% ungroup %>% select(sex) %>% distinct() %>% pull

  lookup <- df %>%
    select(country, cap_adj) %>%
    distinct

  projected <- raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::spread_draws(
      mu_ct[country, year],
      late[country],
      vlate[country]
    ) %>%
    ungroup %>%
    mutate(year = year + baseyear) %>%
    select(-.iteration, -.chain) %>%
    rename(iteration = .draw) %>%
    group_by(country, iteration) %>%
    arrange(year) %>%
    mutate(
      mu3.ct = mu_ct - 2*late,
      mu4.ct = mu_ct - late,
      mu5.ct = mu_ct,
      mu6.ct = mu_ct + vlate,
      mu7.ct = mu_ct + vlate*2,
      mu8.ct = mu_ct + vlate*3
    ) %>%
    select(-mu_ct, -late, -vlate) %>%
    mutate_at(.vars = vars(mu3.ct, mu4.ct, mu5.ct, mu6.ct, mu7.ct, mu8.ct), pnorm) %>%
    left_join(lookup, by="country") %>%
    mutate(mu3.ct = pmax(pmin(mu3.ct + cap_adj, 1), 0),
           mu4.ct = pmax(pmin(mu4.ct + cap_adj, 1), 0),
           mu5.ct = pmax(pmin(mu5.ct + cap_adj, 1), 0),
           mu6.ct = pmax(pmin(mu6.ct + cap_adj, 1), 0),
           mu7.ct = pmax(pmin(mu7.ct + cap_adj, 1), 0),
           mu8.ct = pmax(pmin(mu8.ct + cap_adj, 1), 0)) %>%
    select(country, year, iteration, a3=mu3.ct, a4=mu4.ct, a5=mu5.ct, a6=mu6.ct, a7=mu7.ct, a8=mu8.ct) %>%
    tidyr::gather(series, value, a3, a4, a5, a6, a7, a8) %>%
    group_by(country, year, series) %>%
    tidybayes::point_interval(value, .width = 0.9) %>%
    select(-.width, -.point, -.interval) %>%
    rename(lower = .lower, upper = .upper) %>%
    ungroup %>%
    mutate_at(vars(value, lower, upper), ~ round(., 3)) %>%
    filter(year >= startyear)

  projected %>% mutate(variable = variable, sex = sex)
}

#' Process ABC joint model output.
#'
#' \code{mdl_jt_process} operates on the raw stanfit object outputted by the ABC
#' joint model and summarises the findings in a digestible format. Specifically,
#' variants of the latent true completion rate are computed and presented in
#' 90\% intervals.
#'
#' @param df Input data frame.
#' @param raw_mcmc Raw ABC model MCMC output.
#' @param props Annual population proportions.
#'
#' @return Data frame with columns for country, year, series, value, lower,
#' upper, variable, and sex. Note that series corresponds to one of:
#' \enumerate{
#'     \item projected3t5 (Official completion rate indicator)
#'     \item projected5 (Upper bracket of the indicator range)
#'     \item projected8 (Ultimate completion)
#' }
#'
#' @export
mdl_jt_process <- function(df, raw_mcmc, props) {
  startyear <- min(df$year)
  baseyear <- startyear - 1
  variable <- df %>% ungroup %>% select(variable) %>% distinct() %>% pull

  lookup <- df %>%
    select(country, cap_adj) %>%
    distinct

  projected <- raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::spread_draws(
      mu_ct[sex, country, year],
      late[country, sex],
      vlate[country, sex]
    ) %>%
    ungroup %>%
    mutate(year = year + baseyear) %>%
    select(-.iteration, -.chain) %>%
    rename(iteration = .draw) %>%
    group_by(country, sex, iteration) %>%
    arrange(year) %>%
    mutate(
      mu5.ct = mu_ct,
      mu4.ct = dplyr::lead(mu_ct, 1) - late,
      mu3.ct = dplyr::lead(mu_ct, 2) - 2*late,
      mu8.ct = mu_ct + vlate * 3
    ) %>%
    select(-mu_ct, -late, -vlate) %>%
    tidyr::pivot_longer(cols = starts_with("mu"), names_to = "indicator", values_to = "value") %>%
    tidyr::pivot_wider(names_from = "sex", values_from = "value") %>%
    left_join(props, by = c("country", "year")) %>%
    mutate(total = female*fprop + male*mprop) %>%
    select(-mprop, -fprop) %>%
    tidyr::pivot_longer(cols = c("female", "male", "total"), names_to = "sex", values_to = "value") %>%
    tidyr::pivot_wider(names_from = "indicator", values_from = "value") %>%
    mutate_at(.vars = vars(mu5.ct, mu4.ct, mu3.ct, mu8.ct), pnorm) %>%
    left_join(lookup, by="country") %>%
    mutate(mu5.ct = pmax(pmin(mu5.ct + cap_adj, 1), 0),
           mu4.ct = pmax(pmin(mu4.ct + cap_adj, 1), 0),
           mu3.ct = pmax(pmin(mu3.ct + cap_adj, 1), 0),
           mu8.ct = pmax(pmin(mu8.ct + cap_adj, 1), 0)) %>%
    select(-cap_adj) %>%
    mutate(projected3t5 = (mu3.ct + mu4.ct + mu5.ct)/3) %>%
    rename(projected5 = mu5.ct, projected8 = mu8.ct) %>%
    select(country, year, sex, iteration, projected3t5, projected5, projected8) %>%
    tidyr::gather(series, value, projected3t5, projected5, projected8) %>%
    group_by(country, year, sex, series) %>%
    tidybayes::point_interval(value, .width = 0.9) %>%
    select(-.width, -.point, -.interval) %>%
    rename(lower = .lower, upper = .upper) %>%
    ungroup %>%
    mutate_at(vars(value, lower, upper), ~ round(., 3)) %>%
    filter(year >= startyear)

  projected %>% mutate(variable = variable)
}

#' ABC model PSIS-LOO.
#'
#' Performs PSIS-LOO on the ABC model output.
#'
#' @param mcmc Raw MCMC output.
#'
#' @return LOO object.
#'
#' @export
loo_output <- function(mcmc) {
  lls <- as.array(mcmc, pars="ll")
  loo::loo(lls, r_eff = loo::relative_eff(exp(lls)), save_psis = TRUE)
}

#' Determine influential points.
#'
#' \code{loo_influential} uses the estimated Pareto-k values from the PSIS-LOO
#' method to determine which points are influential to model estimation. Note that
#' due to the model structure, some observations will by design produce high Pareto-k
#' values such as a3 points subject to late completion. While PSIS-LOO may not generate
#' a usable LOO estimate, it does point towards interesting observations.
#'
#' @param df Input data frame.
#' @param loo1 LOO object.
#' @param threshold Pareto-K influential threshold
#'
#' @return The set of observations determined to be highly influential to parameter
#' estimation.
#' @export
loo_influential <- function(df, loo1, threshold = 0.7) {
  df %>%
    bind_cols(data.frame(Pareto.k = loo1$diagnostics$pareto_k)) %>%
    filter(Pareto.k > threshold) %>%
    arrange(-Pareto.k)
}

#' Plot ABC LOO-PIT
#'
#' Plots the LOO-PIT for the ABC model output
#'
#' @param df Input data frame.
#' @param raw_mcmc ABC model output.
#' @param loo1 LOO object.
#'
#' @export
loo_pit <- function(df, raw_mcmc, loo1) {
  yrepl <- rstan::extract(raw_mcmc)[["yrepl"]]
  bayesplot::ppc_loo_pit_overlay(yrep = yrepl, y = qnorm(pmax(pnorm(-2.5), pmin(pnorm(2.5), df$value))),
                                 lw = stats::weights(loo1$psis_object))
}


#' Plot ABC LOO-PIT
#'
#' Plots the LOO-PIT QQ for the ABC model output
#'
#' @param df Input data frame.
#' @param raw_mcmc ABC model output.
#' @param loo1 LOO object.
#'
#' @export
loo_qq <- function(df, raw_mcmc, loo1) {
  yrepl <- rstan::extract(raw_mcmc)[["yrepl"]]
  bayesplot::ppc_loo_pit_qq(yrep = yrepl, y = qnorm(pmax(pnorm(-2.5), pmin(pnorm(2.5), df$value))),
                            lw = stats::weights(loo1$psis_object), size = 1)
}
