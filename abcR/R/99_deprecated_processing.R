#' Extract observed values (Deprecated).
#'
#' @param df Input data frame.
#' @param cutoff Extreme value threshold.
mdl_observed <- function(df, cutoff=0.98) {
  input <- df %>%
    ungroup %>%
    mutate(raw_value = value,
           was_capped = !between(raw_value, 1-cutoff, cutoff),
           value = pmax(1-cap, pmin(cap, value)),
           capped_slice = ifelse(was_capped, raw_value - value, NA),
           value = qnorm(value))

  input %>% select(country, year, survey, raw_value)
}

# readd(cr_input_prim_total) %>%
#   ungroup %>%
#   mutate(raw_value = value,
#          was_capped = !between(raw_value, 0.02, 0.98),
#          value = pmax(0.02, pmin(0.98, value)),
#          capped_slice = ifelse(was_capped, raw_value - value, NA), value = qnorm(value)) %>%
#   group_by(survey) %>%
#   mutate(survey.cap = max(was_capped)) %>%
#   ungroup %>%
#   group_by(country) %>% read
#   mutate(n_survey = n_distinct(survey), country.cap = max(was_capped), any.uncap = min(survey.cap)) %>%
#   mutate(cap.case = case_when(country.cap == 0 ~ 1,
#                               n_survey == 1 & country.cap == 1 ~ 2,
#                               n_survey > 1 & any.uncap == 1 ~ 3,
#                               TRUE ~ 4)) %>%
#   ungroup


#' Process ABC model output.
#'
#' @param df Input data frame.
#' @param raw_mcmc Raw MCMC output.
#' @param rounds Survey rounds data frame.
#'
#' @export
mdl_process2 <- function(df, raw_mcmc, rounds) {
  baseyear    <- min(df$year) - 1
  inputs <-   df %>%
    ungroup %>%
    group_by(country, survey) %>%
    mutate(se_temp = ifelse(is.finite(se_q), se_q, NA),
           se_mean = mean(se_temp, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      raw_value = value,
      was_capped = !between(raw_value, 0.02, 0.98),
      value = pmax(0.02, pmin(0.98, value)),
      capped_slice = ifelse(was_capped, raw_value - value, NA),
      value = qnorm(value),
      se    = case_when(se_q == Inf & !is.na(se_mean) ~ se_mean,
                        se_q == Inf & is.na(se_mean) ~ 0.35,
                        TRUE        ~ se_q)
      # ,index = 1:n()
    )
  start_yrs   <- summarise(group_by(df, country), start = min(year)) %>% mutate(start = baseyear + 1)

  pars <-
    raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::gather_draws(
      mu_ct[country, year],
      beta_s[survey],
      #beta_s_scaled[survey],
      drift[country],
      late[country],
      mult5err[country],
      #sigma,
      #sigma_beta_s,
      sigma_s[survey],
      vlate[country]
    ) %>%
    ungroup %>%
    mutate(year = year + baseyear) %>%
    split(.$.variable) %>%
    purrr::map(~ {select(.x, -.iteration, -.chain) %>%
        rename(iteration = .draw, variable = .variable, value = .value) %>%
        select_if(function(x){!all(is.na(x))})})

  pars_wide_prj <-
    raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::spread_draws(
      mu_ct[country, year],
      #drift[country],
      late[country],
      mult5err[country],
      vlate[country]
    ) %>%
    ungroup %>%
    mutate(year = year + baseyear) %>%
    select(-.iteration, -.chain) %>%
    rename(iteration = .draw)

  pars_wide_fit <-
    raw_mcmc %>%
    tidybayes::recover_types(df) %>%
    tidybayes::spread_draws(
      beta_s[survey],
      #beta_s_scaled[survey],
      #sigma,
      #sigma_beta_s,
      sigma_s[survey]
    ) %>%
    ungroup %>%
    select(-.iteration, -.chain) %>%
    rename(iteration = .draw)

  fit <-
    inputs %>%
    left_join(pars_wide_prj) %>%
    # left_join(inputs) %>%
    left_join(pars_wide_fit) %>%
    mutate(
      yhat =    mu_ct
      + beta_s
      - late * obsage/2
      + vlate * pmin(3, recondist)
      - mult5err * truage5mlt,
      sd = sqrt(se^2 + sigma_s^2),
      yrepl = rnorm(dplyr::n(), mean = yhat, sd = sd),
      ll = dnorm(log = TRUE, x = value, mean = yhat, sd = sd)
    ) %>%
    select(iteration, index, country, year, survey, survey_series, variable, value, age, sex,
           truage5mlt, obsage, recondist,
           yhat, sd, yrepl, ll)

  projected <-
    pars_wide_prj %>%
    group_by(country, iteration) %>%
    arrange(year) %>%
    mutate(
      mu5.ct = mu_ct,
      mu4.ct = dplyr::lead(mu_ct, 1) - late/2,
      mu3.ct = dplyr::lead(mu_ct, 2) - late,
      mu8.ct = mu_ct + vlate * 3
    ) %>%
    select(-mu_ct, -late, -vlate, -mult5err) %>%
    mutate_at(.vars = vars(mu5.ct, mu4.ct, mu3.ct, mu8.ct), pnorm) %>%
    mutate(
      projected3t5 = (mu3.ct + mu4.ct + mu5.ct)/3,
      projected5   = mu5.ct,
      projected8   = mu8.ct
    ) %>%
    select(-mu5.ct, -mu4.ct, -mu3.ct, -mu8.ct) %>%
    tidyr::gather(variable, value, projected3t5, projected5, projected8) %>%
    group_by(country, year, variable) %>%
    tidybayes::point_interval(value, .width = 0.9) %>%
    select(-.width, -.point, -.interval) %>%
    rename(lower = .lower, upper = .upper) %>%
    ungroup %>%
    mutate_at(vars(value, lower, upper), ~ round(., 3)) %>%
    left_join(start_yrs, by = 'country') %>%
    filter(year >= start) %>%
    select(-start)

  forLOO <-
    list(obs = inputs$value,
         ll =
           select(fit, ll, iteration, index) %>%
           na.omit %>%
           arrange(iteration, index) %>%
           {t(matrix(data = .$ll, ncol = max(.$iteration)))},
         reps =
           select(fit, yrepl, iteration, index) %>%
           na.omit %>%
           arrange(iteration, index) %>%
           {t(matrix(data = .$yrepl, ncol = max(.$iteration)))}
    )

  #pars_summary <-
  #  c(pars,
  #    list(beta_s_type =
  #    {pars$beta_s_scaled %>%
  #        left_join(mutate(rounds, survey = paste0(country, survey, year))) %>%
  #        select(iteration, survey, variable, value, type = round) %>%
  #        group_by(iteration, type) %>%
  #        summarize(value = mean(value, na.rm = TRUE))})) %>%
  #  map(~ {group_by_at(.x, .vars = vars(-iteration, -value)) %>%
  #      tidybayes::point_interval(value, .width = 0.9) %>%
  #      select(-.width, -.point, -.interval) %>%
  #      rename(lower = .lower, upper = .upper) %>%
  #      mutate_at(vars(value, lower, upper), ~ round(., 3))})

  list(crs = c(list(observed     = select(df,  country, year, survey, value)
                    , fitted     = select(fit, index, country, year, survey, value = yhat)
                    , replicated = select(fit, index, country, year, survey, value = yrepl)
  ),
  split(select(projected, -variable), projected$variable)),
  #pars = pars_summary,
  fit = fit,
  logliks = forLOO#,
  #pars_iter = pars[c('drift')]
  )
}
