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
  } else if (type == "ran") {
    df %>%
      anti_join(df %>%
                  filter(year != min(year)) %>%
                  group_by(survey) %>%
                  mutate(n = n()) %>%
                  filter(n > 5) %>%
                  select(-n) %>%
                  slice_sample(n=2))
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

#' Model-based ABC test.
#'
#' \code{test_abc2} uses the results of the ABC model executed using the input
#' set to estimate the observations described by the test set. The differences
#' are summarised by the mean squared error, mean absolute deviation, and
#' median absolute deviation. Unlike \code{test_abc}, this function uses a
#' Stan model's generated quantities block to execute tests instead of
#' relying on posterior sample manipulation.
#'
#' @param input Input set data frame.
#' @param tests Test set data frame.
#' @param raw_mcmc Raw Stan model output.
#' @param model Generated quantities model.
#' @param type Test type.
#'
#' @family testing functions
#'
#' @return Returns a summary table with the model performance.
#'
#' @export
test_abc2 <- function(input, tests, raw_mcmc, model, type) {
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

  print(model)

  data <- bind_rows(input %>% mutate(test = 0),
                    tests %>% mutate(test = 1)) %>%
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
  data[["iters"]] <- dim(data[["mu_ct"]])[1]
  data[["late"]] <- get_parsamps(raw_mcmc, "late") %>% replace_na(0)
  data[["vlate"]] <- get_parsamps(raw_mcmc, "vlate") %>% replace_na(0)
  data[["cibounds"]] <- floor(data[["iters"]]*c(0.025,0.05,0.1,0.9,0.95,0.975))
  data[["mult5err"]] <- get_parsamps(raw_mcmc, "mult5err")
  if (is.null(data[["mult5err"]])) {
    data[["mult5err"]] <- matrix(data=0, nrow=data[["iters"]], ncol=data[["n_country"]])
  }

  if (type == "lso") {
    if (is.null(get_parsamps(raw_mcmc, "sigma_s"))) {
      data[["sigma_s"]] <- get_parsamps(raw_mcmc, "sigma_c")
    } else {
      data[["sigma_s"]] <- get_nsvar(input, raw_mcmc)
    }
  } else {
    data[["sigma_s"]] <- get_parsamps(raw_mcmc, "sigma_s")
    if (is.null(data[["sigma_s"]])) {
      temp1 <- get_parsamps(raw_mcmc, "sigma_c")
      temp2 <- bind_rows(input %>% mutate(test = 0),
                         testss %>% mutate(test = 1)) %>%
        ungroup %>%
        select(country, survey) %>%
        distinct %>%
        arrange(survey) %>%
        mutate(country = as.numeric(as.factor(country))) %>%
        select(country) %>%
        pull

      temp3 <- matrix(data=0, nrow=data[["iters"]], ncol=data[["n_survey"]])
      for (i in 1:data[["n_survey"]]) {
        temp3[, i] <- temp1[, temp2[i]]
      }
      data[["sigma_s"]] <- temp3
    }

    data[["beta_s"]] <- get_parsamps(raw_mcmc, "beta_s")
    if (is.null(data[["beta_s"]])) {
      data[["beta_s"]] <- matrix(data=0, nrow=data[["iters"]], ncol=data[["n_survey"]])
    }

  }

  mod <- stan(model, seed = 2020, chains = 1, data = data,
              warmup = 0, iter = 10, algorithm = "Fixed_param")

  bind_cols(
    tsmall,
    rstan::summary(mod, pars=c("yhat"), probs=c("0.5"))$summary %>% as.data.frame
  ) %>%
    left_join(lookup, by = c("country")) %>%
    mutate(pred = pmax(pmin(pnorm(`50%`) + cap_adj, 1), 0)) %>%
    mutate(diff = pred - value, diff_sq = diff^2) %>%
    summarise(mse.100 = round(100 * mean(diff_sq, na.rm = TRUE), 3),
              mae.100 = round(100 * mean(abs(diff), na.rm = TRUE), 3),
              mbe.100 = round(100 * mean(diff, na.rm = TRUE), 3),
              mmad.100 = round(100 * median(abs(diff), na.rm = TRUE), 3)) %>%
    mutate(model = "ABC", variable = variable, sex = sex) %>%
    select(model, level = variable, sex, everything())

  # mcmc %>%
  #   tidybayes::recover_types(input) %>%
  #   tidybayes::spread_draws(
  #     mu_ct[country, year],
  #     late[country],
  #     vlate[country],
  #     mult5err[country]
  #   ) %>%
  #   ungroup %>%
  #   mutate(year = year + baseyear) %>%
  #   select(-.iteration, -.chain) %>%
  #   rename(iteration = .draw) %>%
  #   inner_join(tsmall, by = c("country", "year")) %>%
  #   left_join(lookup, by = c("country")) %>%
  #   mutate(pred = mu_ct + vlate * recondist3 - late * halfobsage - mult5err * truage5mlt,
  #          pred = pmax(pmin(pnorm(pred) + cap_adj, 1), 0)) %>%
  #   select(country, year, iteration, value, pred, recondist) %>%
  #   group_by(country, year, value, recondist) %>%
  #   tidybayes::point_interval(pred, .width = 0.9) %>%
  #   ungroup %>%
  #   select(-.width, -.point, -.interval) %>%
  #   rename(lower = .lower, upper = .upper) %>%
  #   filter(recondist == min(recondist)) %>%
  #   group_by(country) %>%
  #   summarise(value = mean(value), pred = mean(pred)) %>%
  #   mutate(diff = pred - value, diff_sq = diff^2) %>%
  #   summarise(mse.100 = round(100 * mean(diff_sq, na.rm = TRUE), 3),
  #             mae.100 = round(100 * mean(abs(diff), na.rm = TRUE), 3),
  #             mbe.100 = round(100 * mean(diff, na.rm = TRUE), 3),
  #             mmad.100 = round(100 * median(abs(diff), na.rm = TRUE), 3)) %>%
  #   mutate(model = "ABC", variable = variable, sex = sex) %>%
  #   select(model, level = variable, sex, everything())
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
    data[["iters"]] <- dim(data[["mu_ct"]])[1]
    data[["late"]] <- get_parsamps(raw_mcmc, "late") %>% replace_na(0)
    data[["vlate"]] <- get_parsamps(raw_mcmc, "vlate") %>% replace_na(0)
    data[["cibounds"]] <- floor(data[["iters"]]*c(0.025,0.05,0.1,0.9,0.95,0.975))
    data[["mult5err"]] <- get_parsamps(raw_mcmc, "mult5err")
    if (is.null(data[["mult5err"]])) {
      data[["mult5err"]] <- matrix(data=0, nrow=data[["iters"]], ncol=data[["n_country"]])
    }

    if (type == "lso") {
      if (is.null(get_parsamps(raw_mcmc, "sigma_s"))) {
        data[["sigma_s"]] <- get_parsamps(raw_mcmc, "sigma_c")
      } else {
        data[["sigma_s"]] <- get_nsvar(input, raw_mcmc)
      }
    } else {
      data[["sigma_s"]] <- get_parsamps(raw_mcmc, "sigma_s")
      if (is.null(data[["sigma_s"]])) {
        temp1 <- get_parsamps(raw_mcmc, "sigma_c")
        temp2 <- bind_rows(input %>% mutate(test = 0),
                           testset %>% mutate(test = 1)) %>%
          ungroup %>%
          select(country, survey) %>%
          distinct %>%
          arrange(survey) %>%
          mutate(country = as.numeric(as.factor(country))) %>%
          select(country) %>%
          pull

        temp3 <- matrix(data=0, nrow=data[["iters"]], ncol=data[["n_survey"]])
        for (i in 1:data[["n_survey"]]) {
          temp3[, i] <- temp1[, temp2[i]]
        }
        data[["sigma_s"]] <- temp3
      }

      data[["beta_s"]] <- get_parsamps(raw_mcmc, "beta_s")
      if (is.null(data[["beta_s"]])) {
        data[["beta_s"]] <- matrix(data=0, nrow=data[["iters"]], ncol=data[["n_survey"]])
      }

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
