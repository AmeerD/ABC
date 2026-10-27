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
#'   Must contain a \code{cap_adj} column; see \code{\link{add_cap_adj}}.
#' @param cutoff Extreme value cap (DEPRECATED, unused; use \code{\link{add_cap_adj}}).
#' @param end_year Last calendar year of the latent series. The model horizon
#'   (\code{nyears}) is \code{end_year - (min(df$year) - 1)}. If \code{NULL}
#'   (the historical behaviour) it is \code{current year + 5}, which makes the
#'   model dimensions depend on the date the code is run; pass it explicitly for
#'   reproducible runs. Note that the \code{projected3t5} series produced by
#'   \code{\link{mdl_process}} ends at \code{end_year - 2}.
#' @param ... Pass through parameters.
#'
#' @return Outputs a list formatted as input into the ABC Stan model.
#'
#' @export
mdl_prep <- function(df, cutoff = 0.98, end_year = NULL, ...) {
  if (!"cap_adj" %in% names(df)) {
    stop("`df` has no `cap_adj` column. Add it with add_cap_adj() ",
         "(or split_abc_groups()) before calling mdl_prep().", call. = FALSE)
  }
  baseyear    <- min(df$year) - 1
  if (is.null(end_year)) {
    end_year <- lubridate::year(Sys.Date()) + 5
    message("mdl_prep(): `end_year` not supplied; using current year + 5 = ",
            end_year, ". Pass `end_year` explicitly for reproducible runs.")
  }

  out <- df %>%
    #
    group_by(country) %>%
    mutate(haslate = ifelse(median(value) >= 0.95, 0, 1)) %>%
    mutate(
      se_q2 = se/dnorm(qnorm(value)),
      value = value - cap_adj
    ) %>%
    ungroup %>%


    #group_by(country, survey) %>%
    #mutate(se_temp = ifelse(is.finite(se_q), se_q, NA),
    #       se_mean = mean(se_temp, na.rm = TRUE)) %>%
    #ungroup() %>%
    #
    mutate(
      #year = year - min(year) + 1,
      value = qnorm(value)#,
      #se     = case_when(se_q == Inf & !is.na(se_mean) ~ se_mean,
      #                   se_q == Inf & is.na(se_mean) ~ 0.35,
      #                   TRUE        ~ se_q)
    ) %>%
    #select(-se_temp, -se_mean) %>%

    tidybayes::compose_data() %>%
    within(., {
      year = year - baseyear
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
      average = data.frame(country, value) %>%
        mutate(value = pnorm(value)) %>%
        group_by(country) %>%
        summarise(avg = mean(value)) %>%
        arrange(country) %>%
        select(avg) %>%
        pull()
      varcat = data.frame(country, value) %>%
        group_by(country) %>%
        summarise(mn = mean(value), md = median(value)) %>%
        mutate(vc = ifelse(mn<=-2.5 | mn >= 2.5 | md<=-2.5 | md >= 2.5, 2, 1)) %>%
        arrange(country) %>%
        select(vc) %>%
        pull()
    })
  out$nyears <- as.integer(end_year - baseyear)
  if (out$nyears < max(out$year)) {
    stop("`end_year` (", end_year, ") is earlier than the last observed ",
         "reconstructed year (", max(out$year) + baseyear, ").", call. = FALSE)
  }
  out
}

#' Additional ABC model inputs for joint modelling.
#'
#' \code{mdl_joint_prep} appends a matrix of population proportions to the pre-processed
#' model input thus completing the input list for the joint model.
#'
#' @param inputs Pre-processed ABC model input.
#' @param jpars Optional joint modelling parameters.
#' @param ... Pass through parameters.
#'
#' @return Outputs a list formatted as input into the ABC joint Stan model.
#'
#' @export
mdl_joint_prep <- function(inputs, jpars = NULL, ...) {
  if (is.null(jpars)) {
    inputs
  } else {
    c(inputs, jpars)
  }
}

#' Compile ABC model
#'
#' \code{mdl_compile} compiles the given stan model.
#'
#' @param model \code{NULL} (default) to compile the ABC independent model shipped
#'   with the package (\code{\link{abc_model_file}}), a path to a \code{.stan}
#'   file, or Stan code as a character vector (e.g. from \code{readLines()}).
#' @param ... Passed to \code{\link[rstan]{stan_model}} (e.g. \code{auto_write}).
#'
#' @return Compiled stan model.
#' @export
mdl_compile <- function(model = NULL, ...) {
  if (is.null(model)) model <- abc_model_file()
  smod <- if (length(model) == 1L && !grepl("\n", model) && file.exists(model)) {
    stan_model(file = model, ...)
  } else {
    stan_model(model_code = paste(model, collapse = "\n"), ...)
  }
  message("Model compiled.")
  smod
}

#' Initialize the drift parameter.
#'
#' Due to the EMG prior on the drift parameter, the standard Unif(-2,2)
#' can produce values with effective probability 0 due to machine
#' precision limitations. Thus, we consider Unif(-0.5,2) instead.
#'
#' @param ncountry Number of countries.
#' @param nyears Number of years.
#' @param joint Boolean indicating if this is input for the joint model.
#'
#' @return Initial values for the drift parameter
init_drift <- function(ncountry, nyears, joint = FALSE) {
  if (joint) {
    list(struc = runif(ncountry, min = 0, max = 0.2),
         drift = matrix(runif(2*ncountry, min = 0, max = 0.2), ncol=2))
  } else {
    list(drift = runif(ncountry, min = 0, max = 2))
  }
}

#' Execute ABC model.
#'
#' \code{mdl_core} calls the Stan sampler to sample from the ABC model.
#' It is normally called through \code{\link{mdl_run}}.
#'
#' Initial values: unless \code{init} is supplied, one list per chain is drawn
#' with \code{init_drift}, using the seed \code{init_seed + chain_nr - 1}, so a
#' run is determined by \code{seed}, \code{init_seed}, \code{chain_nr} and the
#' data. Runs that fit one chain per call with different \code{chain_nr} (as in
#' \code{scripts/sample_plans.R}) therefore still start from different values.
#' The global random number generator state is left unchanged. (Before 0.2.0 the
#' inits came from the unseeded global RNG, and a single init list was passed
#' to the sampler, so \code{nchains > 1} in one call failed.)
#'
#' @param data Stan data list from \code{\link{mdl_prep}}.
#' @param chain_nr First chain id.
#' @param init \code{"random"} (draw drift inits, see Details), or a list of
#'   per-chain init lists, or anything accepted by \code{\link[rstan]{sampling}}.
#' @param model Compiled ABC Stan model (\code{\link{mdl_compile}}).
#' @param nchains Number of chains.
#' @param nburn Number of warm-up iterations.
#' @param niter Number of post-warm-up iterations.
#' @param nthin Thinning interval.
#' @param seed Sampler seed (historically fixed at 2019).
#' @param init_seed Seed for drawing initial values (defaults to \code{seed});
#'   \code{NULL} draws them from the current RNG state.
#' @param cores Number of chains to run in parallel.
#' @param adapt_delta,max_treedepth NUTS control settings.
#' @param refresh Progress printing interval passed to the sampler.
#' @param ... Ignored (pass-through from \code{\link{mdl_run}}).
#'
#' @return Returns a stanfit object. The initial values used are stored in
#'   \code{attr(fit, "abc_init")}.
mdl_core <- function(data, chain_nr = 1, init = 'random', model, nchains = 1,
                     nburn, niter, nthin = 1, seed = 2019, init_seed = seed,
                     cores = 1, adapt_delta = 0.9, max_treedepth = 12,
                     refresh = max(floor((nburn + niter) / 10), 1), ...) {
  if (missing(model)) stop("Supply a compiled model: model = mdl_compile().", call. = FALSE)
  if (identical(init, 'random')) {
    init <- mdl_draw_inits(data, nchains = nchains,
                           seed = if (is.null(init_seed)) NULL else init_seed + chain_nr - 1)
  }
  fit <- sampling(
    model,
    data = data,
    seed = seed,
    chains = nchains,
    chain_id = seq.int(chain_nr, length.out = nchains),
    cores = cores,
    init = init,
    warmup = nburn,
    iter   = nburn + niter,
    thin   = nthin,
    refresh = refresh,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )
  attr(fit, "abc_init") <- init
  fit
}

#' Draw initial values for the ABC model.
#'
#' @param data Stan data list from \code{\link{mdl_prep}}.
#' @param nchains Number of chains.
#' @param seed Seed for the random draws, or \code{NULL} to use the current RNG
#'   state. With a seed, the global RNG state is restored afterwards.
#'
#' @return A list of \code{nchains} init lists.
#' @export
mdl_draw_inits <- function(data, nchains = 1, seed = NULL) {
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
            else rm(".Random.seed", envir = globalenv()))
    set.seed(seed)
  }
  lapply(seq_len(nchains), function(i)
    init_drift(data$n_country, data$nyears, !is.null(data$mprop)))
}

#' Run the ABC model.
#'
#' \code{mdl_run} pre-processes the input data using \code{\link{mdl_prep}}
#' and subsequently samples from the ABC model using \code{\link{mdl_core}}.
#'
#' @param df Input data frame (one level/sex group, with \code{cap_adj}).
#' @param ... Passed to \code{\link{mdl_prep}} (e.g. \code{end_year}) and
#'   \code{\link{mdl_core}} (e.g. \code{model}, \code{nchains}, \code{nburn},
#'   \code{niter}, \code{nthin}, \code{seed}, \code{cores}).
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
