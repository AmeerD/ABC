# Workflow helpers (abcR 0.2.0). They package steps that scripts/ABCtest.R
# performs inline (cap adjustment, splitting into level/sex groups) together
# with input checks, multi-chain runs and convergence diagnostics.

#' Path to the Stan model shipped with abcR.
#'
#' @param name Model name. \code{"ABC_indep"} is the independent ABC model of
#'   Dharamshi et al. (2022) with array declarations updated to the syntax
#'   required by Stan >= 2.33. The statistical model is unchanged.
#'
#' @return Path to the \code{.stan} file.
#' @export
abc_model_file <- function(name = "ABC_indep") {
  path <- system.file("stan", paste0(name, ".stan"), package = "abcR")
  if (!nzchar(path)) stop("No packaged Stan model called ", name, call. = FALSE)
  path
}

#' Example (anonymised) ABC input data.
#'
#' Loads the public test data set distributed with the upstream ABC repository.
#' Country codes are scrambled and values perturbed, so results cannot be
#' compared with published UIS figures.
#'
#' @return The \code{test.dat} data frame.
#' @export
abc_example_data <- function() {
  env <- new.env(parent = emptyenv())
  load(system.file("extdata", "ABCtest.rda", package = "abcR", mustWork = TRUE),
       envir = env)
  as.data.frame(env$test.dat)
}

#' Add the cap adjustment used before probit transformation.
#'
#' Observed completion rates close to 0 or 1 have unstable probit transforms.
#' For each country/level/sex group, if the maximum exceeds \code{cap} all
#' values are shifted down by \code{max - cap}; otherwise, if the minimum is
#' below \code{1 - cap}, they are shifted up by \code{(1 - cap) - min}. The shift
#' is stored in \code{cap_adj}, subtracted in \code{\link{mdl_prep}} and added
#' back in \code{\link{mdl_process}}. This is the rule used in the upstream
#' example script.
#'
#' @param df Input data with \code{country}, \code{variable}, \code{sex}, \code{value}.
#' @param cap Cap (default 0.98).
#' @param both_extremes What to do when one group has values both above
#'   \code{cap} and below \code{1 - cap}; a single shift cannot fix both, so the
#'   default is to stop.
#'
#' @return \code{df} with a \code{cap_adj} column (existing one replaced).
#' @export
add_cap_adj <- function(df, cap = 0.98, both_extremes = c("error", "warn")) {
  both_extremes <- match.arg(both_extremes)
  caps <- df %>%
    ungroup() %>%
    group_by(country, variable, sex) %>%
    summarise(upper = max(value), lower = min(value), .groups = "drop")
  bad <- caps %>% filter(upper > cap & lower < 1 - cap)
  if (nrow(bad)) {
    msg <- paste0("Values both above ", cap, " and below ", 1 - cap,
                  " within: ", paste(bad$country, bad$variable, bad$sex,
                                     collapse = "; "))
    if (both_extremes == "error") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  caps <- caps %>%
    mutate(cap_adj = case_when(upper > cap ~ upper - cap,
                               lower < (1 - cap) ~ lower - (1 - cap),
                               TRUE ~ 0)) %>%
    select(country, variable, sex, cap_adj)
  df %>%
    select(-any_of("cap_adj")) %>%
    left_join(caps, by = c("country", "variable", "sex"))
}

#' Split ABC input into the nine level/sex model groups.
#'
#' @param df Input data (e.g. \code{\link{abc_example_data}()}).
#' @param cap Passed to \code{\link{add_cap_adj}}.
#' @param ... Passed to \code{\link{add_cap_adj}}.
#'
#' @return Named list with elements among \code{pt, pm, pf, lt, lm, lf, ut, um,
#'   uf} (p/l/u = primary/lower/upper secondary; t/m/f = total/male/female).
#'   Groups without data are omitted.
#' @export
split_abc_groups <- function(df, cap = 0.98, ...) {
  df <- as.data.frame(ungroup(df))
  check_abc_input(df)
  df <- add_cap_adj(df, cap = cap, ...)
  spec <- expand.grid(sex = c("total", "male", "female"),
                      variable = c("prim", "lsec", "usec"),
                      stringsAsFactors = FALSE)
  spec$group <- paste0(substr(spec$variable, 1, 1), substr(spec$sex, 1, 1))
  out <- lapply(seq_len(nrow(spec)), function(i)
    df[df$variable == spec$variable[i] & df$sex == spec$sex[i], , drop = FALSE])
  names(out) <- spec$group
  out[vapply(out, nrow, integer(1)) > 0]
}

#' Check raw ABC input data.
#'
#' Stops with an informative message when fields are missing, non-finite or
#' outside the domain required by the probit transformation.
#'
#' @param df Input data.
#' @return \code{df}, invisibly.
#' @export
check_abc_input <- function(df) {
  needed <- c("country", "variable", "sex", "survey", "year", "value", "se",
              "obsage", "recondist", "truage5mlt")
  miss <- setdiff(needed, names(df))
  if (length(miss)) stop("Missing fields: ", paste(miss, collapse = ", "), call. = FALSE)
  if (anyNA(df[needed])) stop("Missing values in required fields.", call. = FALSE)
  for (f in c("year", "value", "se", "obsage", "recondist", "truage5mlt"))
    if (!is.numeric(df[[f]]) || any(!is.finite(df[[f]])))
      stop("Field `", f, "` must be finite numeric.", call. = FALSE)
  if (any(df$value <= 0 | df$value >= 1))
    stop("`value` must lie strictly inside (0, 1).", call. = FALSE)
  if (any(df$se <= 0)) stop("`se` must be positive.", call. = FALSE)
  if (any(df$year != floor(df$year))) stop("`year` must be whole numbers.", call. = FALSE)
  if (!all(df$variable %in% c("prim", "lsec", "usec")))
    stop("`variable` must be one of prim, lsec, usec.", call. = FALSE)
  if (!all(df$sex %in% c("total", "male", "female")))
    stop("`sex` must be one of total, male, female.", call. = FALSE)
  invisible(df)
}

#' Check a prepared Stan data list.
#'
#' The ABC Stan model assumes that survey ids form consecutive blocks within
#' country (it builds sum-to-zero survey effects from \code{survey_count}).
#' This function verifies that assumption and the transformed values.
#'
#' @param dat Output of \code{\link{mdl_prep}}.
#' @return \code{dat}, invisibly.
#' @export
mdl_check_data <- function(dat) {
  if (any(!is.finite(dat$value)) || any(!is.finite(dat$se_q2)) || any(dat$se_q2 <= 0))
    stop("Non-finite transformed values or standard errors.", call. = FALSE)
  if (dat$nyears < max(dat$year)) stop("Latent horizon shorter than data.", call. = FALSE)
  if (sum(dat$survey_count) != dat$n_survey)
    stop("survey_count does not sum to n_survey.", call. = FALSE)
  if (!setequal(dat$country, seq_len(dat$n_country)) ||
      !setequal(dat$survey, seq_len(dat$n_survey)))
    stop("Country or survey ids are not contiguous.", call. = FALSE)
  ends <- cumsum(dat$survey_count)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  for (i in seq_len(dat$n_country)) {
    if (!setequal(unique(dat$survey[dat$country == i]), seq.int(starts[i], ends[i])))
      stop("Survey ids are not arranged in country blocks (country ", i, ").",
           call. = FALSE)
  }
  invisible(dat)
}

#' Parameters that are constant by construction.
#'
#' Countries without (very) late completion get \code{late}/\code{vlate} fixed
#' at 0 and the padding elements of \code{late_short2}/\code{vlate_short2} are 0;
#' the survey effect of a single-survey country is 0 by the sum-to-zero
#' constraint. These have zero posterior variance and undefined Rhat/ESS.
#'
#' @param dat Output of \code{\link{mdl_prep}}.
#' @return Character vector of parameter names.
#' @export
mdl_expected_constants <- function(dat) {
  c(sprintf("late_short2[%d]", sum(dat$haslate2) + 1L),
    sprintf("vlate_short2[%d]", sum(dat$hasvlate2) + 1L),
    sprintf("late[%d]", which(dat$haslate2 == 0)),
    sprintf("vlate[%d]", which(dat$hasvlate2 == 0)),
    sprintf("beta_s[%d]", cumsum(dat$survey_count)[dat$survey_count == 1]))
}

#' Convergence and sampler diagnostics for an ABC fit.
#'
#' Computes rank-normalised Rhat and bulk/tail ESS (\pkg{posterior}) for every
#' quantity, counts divergences and tree-depth saturation, and checks that the
#' quantities with zero variance are exactly those expected by construction.
#'
#' @param fit A stanfit from \code{\link{mdl_run}}/\code{\link{mdl_core}}.
#' @param dat The Stan data list used for the fit.
#' @param rhat_limit,ess_limit Thresholds for flagging.
#'
#' @return A list with \code{summary} (one row), \code{draws} (per-quantity
#'   diagnostics), \code{constants} and \code{flagged}.
#' @export
mdl_diagnose <- function(fit, dat, rhat_limit = 1.01, ess_limit = 400) {
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required.", call. = FALSE)
  ds <- posterior::summarise_draws(
    posterior::as_draws_array(as.array(fit)),
    "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")
  is_const <- is.finite(ds$sd) & ds$sd == 0
  constants <- ds[is_const, ]
  unexplained <- constants[!constants$variable %in% mdl_expected_constants(dat) |
                             constants$mean != 0, ]
  rest <- ds[!is_const, ]
  flagged <- rest[!is.finite(rest$rhat) | rest$rhat >= rhat_limit |
                    !is.finite(rest$ess_bulk) | rest$ess_bulk < ess_limit |
                    !is.finite(rest$ess_tail) | rest$ess_tail < ess_limit, ]
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  max_td <- fit@stan_args[[1]]$control$max_treedepth
  if (is.null(max_td)) max_td <- 10
  energy <- vapply(sp, function(x) {
    e <- x[, "energy__"]; mean(diff(e)^2) / stats::var(e)
  }, numeric(1))
  summary <- data.frame(
    chains = length(sp),
    draws = sum(vapply(sp, nrow, integer(1))),
    divergences = sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1))),
    treedepth_hits = sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= max_td), numeric(1))),
    max_rhat = max(rest$rhat, na.rm = TRUE),
    min_ess_bulk = min(rest$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(rest$ess_tail, na.rm = TRUE),
    constant_quantities = nrow(constants),
    unexplained_constants = nrow(unexplained),
    flagged_quantities = nrow(flagged),
    min_energy_ratio = min(energy),
    thin = fit@stan_args[[1]]$thin
  )
  summary$pass <- summary$divergences == 0 && summary$treedepth_hits == 0 &&
    summary$flagged_quantities == 0 && summary$unexplained_constants == 0 &&
    summary$chains > 1
  list(summary = summary, draws = ds, constants = constants, flagged = flagged)
}

#' Run, process and diagnose all ABC groups.
#'
#' A script-free replacement for the upstream example: prepares each group,
#' samples with explicit seeds and horizon, post-processes with
#' \code{\link{mdl_process}} and computes \code{\link{mdl_diagnose}}. When
#' \code{output_dir} is given each group's fit, Stan data, inits and settings are
#' saved, and an interrupted run resumes; a saved fit is reused only if its
#' stored Stan data and settings are identical to the current ones.
#'
#' @param inputs Named list from \code{\link{split_abc_groups}}.
#' @param model Compiled model; compiled from \code{\link{abc_model_file}} if NULL.
#' @param end_year Last latent year (see \code{\link{mdl_prep}}).
#' @param chains,warmup,iter,thin Sampling settings (\code{iter} counts
#'   post-warm-up iterations, as \code{niter} in \code{\link{mdl_core}}).
#' @param seed Base seed; group \code{i} uses \code{seed + i}.
#' @param cores Parallel chains.
#' @param adapt_delta,max_treedepth NUTS control.
#' @param output_dir Optional directory for saved fits.
#' @param groups Subset of \code{names(inputs)} to run.
#'
#' @return A list with \code{output} (combined \code{mdl_process} output),
#'   \code{diagnostics} (one row per group) and \code{fits} (file paths or
#'   stanfit objects).
#' @export
abc_run_groups <- function(inputs, model = NULL, end_year, chains = 4, warmup = 5000,
                           iter = 5000, thin = 5, seed = 2019, cores = 1,
                           adapt_delta = 0.9, max_treedepth = 12,
                           output_dir = NULL, groups = names(inputs)) {
  if (missing(end_year)) stop("Supply `end_year` explicitly.", call. = FALSE)
  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  settings <- list(chains = chains, warmup = warmup, iter = iter, thin = thin,
                   adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                   end_year = end_year)
  outputs <- list(); diags <- list(); fits <- list()
  for (g in groups) {
    i <- match(g, names(inputs))
    df <- inputs[[g]]
    dat <- mdl_check_data(mdl_prep(df, end_year = end_year))
    group_seed <- seed + i
    fit <- NULL
    if (!is.null(output_dir)) {
      fit_path <- file.path(output_dir, paste0(g, "_fit.rds"))
      meta_path <- file.path(output_dir, paste0(g, "_meta.rds"))
      if (file.exists(fit_path) && file.exists(meta_path)) {
        meta <- readRDS(meta_path)
        if (identical(meta$data, dat) && identical(meta$settings, settings) &&
            identical(meta$seed, group_seed)) {
          message("Reusing saved fit for ", g)
          fit <- readRDS(fit_path)
        } else {
          stop("Saved fit for ", g, " in ", output_dir, " was made with different ",
               "data or settings. Use another output_dir.", call. = FALSE)
        }
      }
    }
    if (is.null(fit)) {
      if (is.null(model)) model <- mdl_compile()
      message("Sampling group ", g, " (", chains, " chains)")
      fit <- mdl_core(dat, model = model, nchains = chains, nburn = warmup,
                      niter = iter, nthin = thin, seed = group_seed,
                      cores = cores, adapt_delta = adapt_delta,
                      max_treedepth = max_treedepth)
      if (!is.null(output_dir)) {
        saveRDS(fit, fit_path)
        saveRDS(list(data = dat, settings = settings, seed = group_seed,
                     init = attr(fit, "abc_init"),
                     abcR_version = as.character(utils::packageVersion("abcR")),
                     rstan_version = as.character(utils::packageVersion("rstan")),
                     created = Sys.time()), meta_path)
      }
    }
    outputs[[g]] <- mdl_process(df, fit)
    diags[[g]] <- cbind(group = g, mdl_diagnose(fit, dat)$summary)
    fits[[g]] <- if (is.null(output_dir)) fit else fit_path
    rm(fit); gc(verbose = FALSE)
  }
  list(output = bind_rows(outputs), diagnostics = bind_rows(diags), fits = fits)
}
