# simple_planeval <- function(p) eval(parse(text = capture.output(plan_to_code(p))),
#                                     envir = .GlobalEnv)

#' Make input data by survey (Deprecated).
#'
#' @param df Raw data.
#'
#' @export
mk_agg_age_by_survey <- function(df) {
  filter(df, between(age, 10, 45)) %>%
    mutate(variable = level) %>%
    wide_aggregate(categories = c('age'), estimate_se = TRUE)
}

#' Aggregate function (Deprecated).
#'
#' @param df Input data frame.
#' @param as Ages?
#' @param ... Pass through parameters.
#'
#' @export
fn_agg <- function(df, as, ...) {
  df %>%
    mutate(
      variable = level,
      age      = age - age_adjustment
    ) %>%
    semi_join(as, by = c('country', 'year', 'level', 'age')) %>%
    wide_aggregate(...)}

#' Country repeater (Deprecated).
#'
#' Function to prepare country list for slicing.
#'
#' @param clist Vector of countries.
#' @param lower Lower bound country.
#' @param upper Upper bound country.
country_rep <- function(clist, lower, upper) {
  ctemp <- clist[clist >= lower & clist <= upper]
  reps <- c()
  for (i in 1:length(ctemp)) {
    reps <- append(reps, ctemp[-i])
  }
  reps
}

# List of parameters for an older version of the model
pars <- c(
  'beta_s'
  ,'beta_s_scaled'
  ,'drift'
  ,'late_c'
  ,'mu_ct'
  ,'mult5err'
  # ,'sigma'
  ,'sigma_beta_s'
  ,'sigma_s'
  ,'vlate_c'
)

#' Generate initial parameter estimates (Deprecated).
#'
#' @param data Input data list.
#' @param chain_nr NUTS chain number.
#' @param modelfile Path to Stan model.
#' @param nchains Number of chains.
#' @param nthin Chain thinning parameter.
#' @param ... Pass through parameters.
#'
#' @export
mdl_inits <- function(data, chain_nr, modelfile, nchains, nthin, ...) {
  mod <- stan(
    seed = 2019,
    chain_id = chain_nr,
    chains = nchains,
    data = data,
    init_r = 2,
    #pars = pars,
    warmup = 5000,
    iter   = 10000,
    thin   = nthin,
    control = list(adapt_delta = 0.75, max_treedepth = 7),
    model_code = modelfile)

  inits <- as.list(summary(mod)$summary[, 'mean'])
  return(inits)
}

#' Run the initialization model.
#'
#' Pre-process and run a faster initial run of the model.
#'
#' @param df Input data frame.
#' @param ... Pass through parameters.
#'
#' @export
mdl_init <- function(df, ...) {
  df %>%
    mdl_prep(...) %>%
    mdl_inits(...)
}

# mdl_core <- function(data, chain_nr, init = 'random') {
#   mod <- rjags::rj(
#     seed = 2019,
#     chain_id = chain_nr,
#     chains = N_chains,
#     data = data,
#     init = init,
#     pars = pars,
#     warmup = N_burnin,
#     iter   = N_burnin + N_iter,
#     thin   = N_thin,
#     control = list(adapt_delta = 0.85, max_treedepth = 10),
#     model_code = modelfile)
#   return(mod)
# }


