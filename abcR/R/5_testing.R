#' Remove latest survey by country.
#'
#' Removes the latest survey from the data for countries with more than one survey. If
#' the user for any reason does not want to filter out any observations (such as while
#' static branching with drake), \code{remove_latest} can simply act as an identity
#' function.
#'
#' @param df Input data frame.
#' @param type Filter type.
#'
#' @family testing functions
#'
#' @export
remove_latest <- function(df, type = "filter") {
  if (type == "full") {
    df
  } else if (type == "rec") {
    remove_recent(df)
  } else {
    df %>%
      group_by(country) %>%
      mutate(nsurvey = n_distinct(survey),
             latest = max(obsyear)) %>%
      filter(obsyear != latest | nsurvey == 1) %>%
      select(-nsurvey, -latest) %>%
      ungroup
  }
}

#' Remove latest surveys by year.
#'
#' Removes all surveys beyond the specified threshold year
#'
#' @param df Input data frame.
#' @param threshold Cutoff year.
#'
#' @family testing functions
#'
#' @export
remove_recent <- function(df, threshold = 2015) {
  df %>%
    group_by(country) %>%
    mutate(sn = n_distinct(survey)) %>%
    filter(year <= threshold | sn == 1) %>%
    select(-sn)
}

#' Observation level ABC test.
#'
#' \code{test_abc} uses the results of the ABC model executed using the input
#' set to estimate the observations described by the test set. The differences
#' are summarised by the mean squared error, mean absolute deviation, and
#' median absolute deviation.
#'
#' @param input Input set data frame.
#' @param tests Test set data frame.
#' @param mcmc Raw ABC MCMC output.
#'
#' @family testing functions
#'
#' @return Returns a summary table with the model performance.
#'
#' @export
test_abc <- function(input, tests, mcmc) {
  if (nrow(tests) == 0) {
    empty <- input %>%
      select(variable, sex) %>%
      distinct %>%
      mutate(model = "ABC",
             mse.100 = NA,
             mae.100 = NA,
             mbe.100 = NA,
             mmad.100 = NA) %>%
      select(model, everything())
    return(empty)
  }

  startyear <- min(input$year)
  baseyear <- startyear - 1
  variable <- input %>% ungroup %>% select(variable) %>% distinct() %>% pull
  sex <- input %>% ungroup %>% select(sex) %>% distinct() %>% pull

  # lookup <- input %>%
  #   select(country, cap_adj) %>%
  #   distinct

  tsmall <- tests %>%
    select(country, year, truage5mlt, obsage, recondist, value, cap_adj) %>%
    mutate(recondist3 = pmin(3, recondist),
           halfobsage = obsage/2) %>%
    select(-obsage, -recondist)

  mcmc %>%
    tidybayes::recover_types(input) %>%
    tidybayes::spread_draws(
      mu_ct[country, year],
      late[country],
      vlate[country],
      mult5err[country]
    ) %>%
    ungroup %>%
    mutate(year = year + baseyear) %>%
    select(-.iteration, -.chain) %>%
    rename(iteration = .draw) %>%
    inner_join(tsmall, by = c("country", "year")) %>%
    #left_join(lookup, by = c("country")) %>%
    mutate(pred = mu_ct + vlate * recondist3 - late * halfobsage - mult5err * truage5mlt,
           pred = pmax(pmin(pnorm(pred) + cap_adj, 1), 0)) %>%
    select(country, year, iteration, value, pred) %>%
    group_by(country, year, value) %>%
    tidybayes::point_interval(pred, .width = 0.9) %>%
    ungroup %>%
    select(-.width, -.point, -.interval) %>%
    rename(lower = .lower, upper = .upper) %>%
    mutate(diff = pred - value, diff_sq = diff^2) %>%
    summarise(mse.100 = round(100 * mean(diff_sq, na.rm = TRUE), 3),
              mae.100 = round(100 * mean(abs(diff), na.rm = TRUE), 3),
              mbe.100 = round(100 * mean(diff, na.rm = TRUE), 3),
              mmad.100 = round(100 * median(abs(diff), na.rm = TRUE), 3)) %>%
    mutate(model = "ABC", variable = variable, sex = sex) %>%
    select(model, level = variable, sex, everything())
}

#' Survey level ABC test.
#'
#' \code{test_abc2} uses the results of the ABC model executed using the input
#' set to estimate the completion rate indicator of the observation year of
#' each survey in the test set. The differences are summarised by the mean
#' squared error, mean absolute deviation, and median absolute deviation.
#'
#' @param input Input set data frame.
#' @param tests Test set data frame.
#' @param mcmc Raw ABC MCMC output.
#'
#' @family testing functions
#'
#' @return Returns a summary table with the model performance.
#'
#' @export
test_abc2 <- function(input, tests, mcmc) {
  if (nrow(tests) == 0) {
    empty <- input %>%
      select(variable, sex) %>%
      distinct %>%
      mutate(model = "ABC",
             mse.100 = NA,
             mae.100 = NA,
             mbe.100 = NA,
             mmad.100 = NA) %>%
      select(model, everything())
    return(empty)
  }

  startyear <- min(input$year)
  baseyear <- startyear - 1
  variable <- input %>% ungroup %>% select(variable) %>% distinct() %>% pull
  sex <- input %>% ungroup %>% select(sex) %>% distinct() %>% pull

  lookup <- input %>%
    select(country, cap_adj) %>%
    distinct

  tsmall <- tests %>%
    select(country, year, truage5mlt, obsage, recondist, value) %>%
    mutate(recondist3 = pmin(3, recondist),
           halfobsage = obsage/2) %>%
    select(-obsage)

  mcmc %>%
    tidybayes::recover_types(input) %>%
    tidybayes::spread_draws(
      mu_ct[country, year],
      late[country],
      vlate[country],
      mult5err[country]
    ) %>%
    ungroup %>%
    mutate(year = year + baseyear) %>%
    select(-.iteration, -.chain) %>%
    rename(iteration = .draw) %>%
    inner_join(tsmall, by = c("country", "year")) %>%
    left_join(lookup, by = c("country")) %>%
    mutate(pred = mu_ct + vlate * recondist3 - late * halfobsage - mult5err * truage5mlt,
           pred = pmax(pmin(pnorm(pred) + cap_adj, 1), 0)) %>%
    select(country, year, iteration, value, pred, recondist) %>%
    group_by(country, year, value, recondist) %>%
    tidybayes::point_interval(pred, .width = 0.9) %>%
    ungroup %>%
    select(-.width, -.point, -.interval) %>%
    rename(lower = .lower, upper = .upper) %>%
    filter(recondist == min(recondist)) %>%
    group_by(country) %>%
    summarise(value = mean(value), pred = mean(pred)) %>%
    mutate(diff = pred - value, diff_sq = diff^2) %>%
    summarise(mse.100 = round(100 * mean(diff_sq, na.rm = TRUE), 3),
              mae.100 = round(100 * mean(abs(diff), na.rm = TRUE), 3),
              mbe.100 = round(100 * mean(diff, na.rm = TRUE), 3),
              mmad.100 = round(100 * median(abs(diff), na.rm = TRUE), 3)) %>%
    mutate(model = "ABC", variable = variable, sex = sex) %>%
    select(model, level = variable, sex, everything())
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

#' Test model on new data
#'
#' \code{mdl_test} uses the results of a model run to test construct the posterior
#' predictive distributions for new data.
#'
#' @param input Input data frame.
#' @param testset Test data frame.
#' @param raw_mcmc Raw Stan model output.
#' @param model Generated quantities model.
#' @param type Test type.
#'
#' @return Returns a stanfit object.
#' @export
mdl_test <- function(input, testset, raw_mcmc, model, type) {
  if (type == "full") {
    return(raw_mcmc)
  } else {
    print(model)

    data <- bind_rows(input %>% mutate(test = 0),
                      testset %>% mutate(test = 1)) %>%
      mdl_prep()

    idx <- data[["test"]] == 1

    #Remove train data
    data[["n"]] <- sum(data[["test"]])
    data[["year"]] <- data[["year"]][idx]
    data[["country"]] <- data[["country"]][idx]
    data[["survey"]] <- data[["survey"]][idx]
    data[["obsage"]] <- data[["obsage"]][idx]
    data[["recondist"]] <- data[["recondist"]][idx]
    data[["truage5mlt"]] <- data[["truage5mlt"]][idx]
    data[["se_q2"]] <- data[["se_q2"]][idx]
    data[["value"]] <- data[["value"]][idx]

    #Append parameter samples
    data[["mu_ct"]] <- get_parsamps(raw_mcmc, "mu_ct")
    data[["mult5err"]] <- get_parsamps(raw_mcmc, "mult5err")
    data[["late"]] <- get_parsamps(raw_mcmc, "late")
    data[["vlate"]] <- get_parsamps(raw_mcmc, "vlate")
    data[["iters"]] <- dim(data[["late"]])[1]
    data[["cibounds"]] <- floor(data[["iters"]]*c(0.025,0.05,0.1,0.9,0.95,0.975))

    if (type == "lso") {
      data[["sigma_s"]] <- get_nsvar(input, raw_mcmc)
    } else {
      data[["sigma_s"]] <- get_parsamps(raw_mcmc, "sigma_s")
      data[["beta_s"]] <- get_parsamps(raw_mcmc, "beta_s")
    }



    mod <- stan(model, seed = 2020, chains = 1, data = data,
                warmup = 0, iter = 10, algorithm = "Fixed_param")
    return(mod)
  }

}

#' Compute coverage probabilities
#'
#' @param input Input data frame.
#' @param raw_mcmc Raw Stan model output.
#' @param type Test type.
#'
#' @return Coverage probabilities.
#' @export
test_cov <- function(input, raw_mcmc, type) {
  if (type == "full") {
    reps <- raw_mcmc %>%
      tidybayes::gather_draws(yrepl[n]) %>%
      select(-.chain, -.iteration) %>%
      tidybayes::point_interval(.value, .width = c(0.8, 0.9, 0.95)) %>%
      select(-.value, -.point, -.interval) %>%
      rename(lower = .lower, upper = .upper) %>%
      pivot_wider(names_from = .width, values_from = c(lower, upper)) %>%
      ungroup
    input %>%
      bind_cols(reps) %>%
      mutate(qval = qnorm(value - cap_adj)) %>%
      mutate(in80 = (qval >= lower_0.8 & qval <= upper_0.8),
             in90 = (qval >= lower_0.9 & qval <= upper_0.9),
             in95 = (qval >= lower_0.95 & qval <= upper_0.95)) %>%
      ungroup %>%
      summarise(cov80 = sum(in80)/n(),
                cov90 = sum(in90)/n(),
                cov95 = sum(in95)/n())
  } else {
    covs <- rstan::get_posterior_mean(raw_mcmc)
    data.frame(coverage = row.names(covs), covs, row.names=NULL) %>%
      filter(coverage != "lp__") %>%
      select(coverage, mean = mean.chain.1) %>%
      pivot_wider(names_from = "coverage", values_from = "mean")
  }
}
