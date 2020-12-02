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
#' upper, vairable, and sex. Note that series corresponds to one of:
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
      mu4.ct = dplyr::lead(mu_ct, 1) - late/2,
      mu3.ct = dplyr::lead(mu_ct, 2) - late,
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
  bayesplot::ppc_loo_pit_overlay(yrep = yrepl, y = qnorm(pmax(0.02, pmin(0.98, df$value))),
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
  bayesplot::ppc_loo_pit_qq(yrep = yrepl, y = qnorm(pmax(0.02, pmin(0.98, df$value))),
                            lw = stats::weights(loo1$psis_object), size = 1)
}
