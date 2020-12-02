#' Construct ABC model input.
#'
#' \code{mk_dataset_4modeling} is the core ABC model data compilation function.
#' It joins a summary of the surveys of interest with the aggregate survey data
#' to produce a complete data frame spanning all levels and sex combinations.
#' Finally, columns detailing the data reconstruction process are appended.
#'
#' @param rounds Survey wave data frame.
#' @param agedat Aggregate age raw data (agg_age).
#' @param genderdat Aggregate age by gender raw data (agg_age_gender).
#' @param ageref Reference ages (refages).
#'
#' @return Raw input data frame.
#'
#' @export
mk_dataset_4modeling <- function(rounds, agedat, genderdat, ageref) {
  bind_rows(
    mutate(agedat, sex = 'total'),
    mutate(genderdat, sex = ifelse(sex == 0, 'female', 'male')),
    tibble()) %>%
    left_join(select(ageref, -year), by = c('country', 'variable' = 'level')) %>%
    filter(refage <= age & age <= refage + 20) %>%
    left_join(rounds) %>%
    rename(survey_series = round) %>%
    # mutate(value = pmax(0.02, pmin(0.98, value))) %>%
    wide_bkprj_4mdl %>%
    whitelist
}


#' Append survey and reconstruction details.
#'
#' @param df Raw input data frame.
wide_bkprj_4mdl <- function(df) {
  df %>%
    mutate(
      survey     = paste0(country, survey, year),
      truage5mlt = as.numeric((age %% 5) == 0),
      obsage     = pmax(0, 2 + refage - age),
      recondist  = pmax(0, age - refage - 2),
      obsyear    = year,
      year       = year + (2 + refage - age)
    ) %>%
    select(-refage)
}

#' Pre-process ABC model input.
#'
#' \code{mdl_prep} is the core data pre-processing function for the ABC model
#' pipeline. It takes a data frame constructed by \code{\link{mk_dataset_4modeling}},
#' adds on additional pieces needed by the model such as indicators for late and
#' very late completion, and transforms the data from the original space to the
#' probit space. Finally, the data is converted to a list suitable to \code{\link[rstan]{stan}}
#' modeling using the \code{\link[tidybayes]{compose_data}} function.
#'
#' @param df Raw input data frame produced by \code{\link{mk_dataset_4modeling}}.
#' @param cutoff Extreme value cap (DEPRECATED).
#' @param ... Pass through parameters.
#'
#' @return Outputs a list formatted as input into the ABC Stan model.
#'
#' @export
mdl_prep <- function(df, cutoff = 0.98, ...) {
  baseyear    <- min(df$year) - 1
  
  df %>%
    #
    group_by(country) %>%
    mutate(haslate = ifelse(median(value) >= 0.95, 0, 1)) %>%
    mutate(
      value = value - cap_adj,
      se_q2 = se/dnorm(qnorm(value))
    ) %>%
    ungroup %>%
    
    
    #group_by(country, survey) %>%
    #mutate(se_temp = ifelse(is.finite(se_q), se_q, NA),
    #       se_mean = mean(se_temp, na.rm = TRUE)) %>%
    #ungroup() %>%
    #
    mutate(
      # year = year - min(year) + 1,
      value = qnorm(value)#,
      #se     = case_when(se_q == Inf & !is.na(se_mean) ~ se_mean,
      #                   se_q == Inf & is.na(se_mean) ~ 0.35,
      #                   TRUE        ~ se_q)
    ) %>%
    #select(-se_temp, -se_mean) %>%
    
    tidybayes::compose_data() %>%
    within(., {
      year = year - baseyear
      nyears = lubridate::year(Sys.Date()) + 5 - baseyear
      survey_count = data.frame(country, survey) %>%
        distinct() %>%
        count(country) %>%
        arrange(country) %>%
        select(n) %>%
        pull()
      haslate2 = data.frame(country, survey, haslate, value, obsage, recondist) %>%
        mutate(value = pnorm(value)) %>%
        group_by(country) %>%
        mutate(ltype = case_when(
          obsage > 0 ~ "late",
          recondist <= 2 ~ "vlate",
          TRUE ~ "irr"
        )) %>%
        group_by(country, survey, haslate, ltype) %>%
        summarise(lmean = mean(value)) %>%
        tidyr::pivot_wider(names_from = ltype, values_from = lmean) %>%
        select(-irr) %>%
        mutate(late_id = case_when(
          haslate == 1 ~ 1,
          is.na(late) | is.na(vlate) ~ 0,
          late < vlate ~ 1,
          TRUE ~ 0)) %>%
        ungroup() %>%
        group_by(country) %>%
        summarise(m = max(late_id)) %>%
        arrange(country) %>%
        select(m) %>%
        pull()
      hasvlate2 = data.frame(country, haslate) %>%
        distinct() %>%
        arrange(country) %>%
        select(haslate) %>%
        pull()
    })
}

#' Additional ABC model inputs for joint modelling.
#'
#' \code{mdl_joint_prep} appends a matrix of population proportions to the pre-processed
#' model input thus completing the input list for the joint model.
#'
#' @param inputs Pre-processed ABC model input.
#' @param props Optional male proportions matrix.
#' @param ... Pass through parameters.
#'
#' @return Outputs a list formatted as input into the ABC joint Stan model.
#'
#' @export
mdl_joint_prep <- function(inputs, props = NULL, ...) {
  if (is.null(props)) {
    inputs
  } else {
    inputs[["mprop"]] = props
    inputs
  }
}

#' Compile ABC model
#'
#' \code{mdl_compile} compiles the given stan model.
#'
#' @param model Stan model file.
#'
#' @return Compiled stan model.
#' @export
mdl_compile <- function(model = modelfile) {
  smod <- stan_model(model_code = model)
  print("Model Compiled!")
  return(smod)
}

#' Initialize the drift parameter.
#'
#' Due to the EMG prior on the drift parameter, the standard Unif(-2,2)
#' can produce values with effective probability 0 due to machine
#' precision limitations. Thus, we consider Unif(-0.5,2) instead.
#'
#' @param ncountry Number of countries.
#'
#' @return Initial values for the drift parameter
init_drift <- function(ncountry, nyears, joint = FALSE) {
  if (joint) {
    list(struc = runif(ncountry, min = -0.2, max = 0.2),
         drift = matrix(runif(2*ncountry, min = -0.2, max = 0.2), nrow=2),
         e_ct = array(runif(2*ncountry*nyears, min=-0.2, max=0.2), dim=c(2,ncountry,nyears)))
  } else {
    list(drift = runif(ncountry, min = -0.5, max = 2))
  }
}

#' Execute ABC model.
#'
#' \code{mdl_core} calls the Stan sampler to sample from the ABC model.
#' It should not be called directly, and instead should only be called
#' through the \code{\link{mdl_run}} function.
#'
#' @param data Input data frame.
#' @param chain_nr NUTS chain number.
#' @param init Initialization method.
#' @param model ABC Stan model.
#' @param ... Pass through parameters.
#'
#' @return Returns a stanfit object.
mdl_core <- function(data, chain_nr, init = 'random', model, nchains, nburn, niter, nthin, ...) {
  init2 <- list(init_drift(data$n_country, data$nyears, !is.null(data$mprop)))
  #smod <- stan_model(model_code = model)
  #print("Model Compiled!")
  mod <- sampling(
    model,
    seed = 2019,
    chain_id = chain_nr,
    chains = nchains,
    data = data,
    init = init2,
    #pars = pars,
    warmup = nburn,
    iter   = nburn + niter,
    thin   = nthin,
    control = list(adapt_delta = 0.8, max_treedepth = 10),
    #model_code = model
  )
  return(mod)
}

#' Run the ABC model.
#'
#' \code{mdl_run} pre-processes the input data using \code{\link{mdl_prep}}
#' and subsequently samples from the ABC model using \code{\link{mdl_core}}.
#'
#' @param df Input data frame.
#' @param ... Pass through parameters.
#'
#' @return Returns a stanfit object.
#'
#' @export
mdl_run <- function(df, ...) {
  df %>%
    mdl_prep(...) %>%
    mdl_joint_prep(...) %>%
    mdl_core(...)
}
