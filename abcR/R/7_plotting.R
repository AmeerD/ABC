# Local lookup table for plotting
label_values <- c(
  "prim" = "Primary",
  "lsec" = "Lower Secondary",
  "usec" = "Upper Secondary",
  "female" = "Female",
  "male" = "Male",
  "total" = "Total"
)

#' Create an ABC input backprojection plot
#'
#' \code{mk_input_plot} plots the input data for a series of countries with
#' separate series for each input survey and identifiers for the age observed
#' and backprojection distance. It is intended to provide a complete picture
#' of the reconstructed surveys for a given country.
#'
#' @param input Input data frame.
#' @param countries Vector of countries to plot.
#'
#' @family plotting functions
#'
#' @export
mk_input_plot <- function(input, countries) {
  input %>%
    filter(country %in% countries) %>%
    ggplot(aes(x=year, y=value)) +
    geom_line(aes(group=survey, colour=survey)) +
    geom_point(aes(shape=as.factor(obsage), alpha = -recondist)) +
    facet_wrap(~country) +
    labs(y = "primary completion") +
    scale_y_continuous(limits = c(0,1), breaks = 0:4 * 0.25) +
    scale_shape_discrete(name = "age observed \n (above nominal graduation age)",
                         labels = c(5, 4, 3)) +
    scale_alpha_continuous(name = "backprojection distance") +
    guides(shape = guide_legend(nrow=1, label.position="bottom", order=1),
           alpha = guide_legend(nrow=1, label.position="bottom", order=2),
           colour = guide_legend(order=3)) +
    theme_minimal()
}

#' Create an ABC input survey plot
#'
#' \code{mk_input_survey_plot} plots the input data for a series of surveys with
#' identifiers for the age observed, backprojection distance, and observations
#' at risk of age misreporting distortions.
#'
#' @param input Input data frame.
#' @param survs Vector of surveys to plot.
#'
#' @family plotting functions
#'
#' @export
mk_input_survey_plot <- function(input, survs) {
  mult5 <- input %>%
    filter(survey %in% survs) %>%
    filter(truage5mlt == 1)

  input %>%
    filter(survey %in% survs) %>%
    ggplot(aes(x=year, y=value)) +
    geom_line(aes(group=survey, colour=survey)) +
    geom_point(data=mult5, aes(x=year, y=value), shape=3, stroke=3, colour="red") +
    geom_point(aes(shape=as.factor(obsage), alpha = -recondist)) +
    facet_wrap(~survey, scales="free_x") +
    labs(y = "primary completion") +
    scale_y_continuous(limits = c(0,1), breaks = 0:4 * 0.25) +
    scale_shape_discrete(name = "age observed \n (above nominal graduation age)",
                         labels = c(5, 4, 3)) +
    scale_alpha_continuous(name = "backprojection distance") +
    guides(shape = guide_legend(nrow=1, label.position="bottom", order=1),
           alpha = guide_legend(nrow=1, label.position="bottom", order=2),
           colour = guide_legend(order=3)) +
    theme_minimal()
}

#' Plot ABC Results by Country
#'
#' \code{plt_country} creates a faceted plot displaying the observed and
#' estimated true completion rates for a given country. Each facet presents
#' a single level/sex combination. Users are given the option to highlight
#' points that were deemed influential in estimating some parameter in the
#' ABC model.
#'
#' @param mdl_proj Projected values data frame.
#' @param mdl_obs Observed values data frame.
#' @param mdl_influential Influential observations data frame.
#' @param use_influential Boolean to indicate whether to plot influential observations.
#' @param target Country code.
#' @param target_yr Final year to plot.
#' @param levels Vector of education levels (prim, lsec, and usec) to plot.
#' @param sexes Vector of sexes (female, male, and total) to plot.
#'
#' @family plotting functions
#'
#' @export
plt_country <- function(mdl_proj, mdl_obs, mdl_influential, target, use_influential = TRUE, target_yr = 2020,
                        levels = c("prim", "lsec", "usec"), sexes = c("female", "male", "total")) {
  output <- mdl_proj %>%
    filter(country == target & year <= target_yr & variable %in% levels & sex %in% sexes) %>%
    mutate(country = countrycode::countrycode(country, 'iso3c', 'country.name')) %>%
    mutate(variable = factor(variable, levels = levels))
  obs <- mdl_obs %>%
    filter(country == target & year <= target_yr & variable %in% levels & sex %in% sexes) %>%
    mutate(country = countrycode::countrycode(country, 'iso3c', 'country.name')) %>%
    mutate(variable = factor(variable, levels = levels))
  influential <- mdl_influential %>%
    filter(country == target & year <= target_yr & variable %in% levels & sex %in% sexes) %>%
    mutate(country = countrycode::countrycode(country, 'iso3c', 'country.name')) %>%
    mutate(variable = factor(variable, levels = levels))
  cht <- output %>% filter(series == "projected8")
  ind <- output %>% filter(series == "projected3t5")

  if (!use_influential) {
    influential <- influential %>% anti_join(influential)
  }

  ggplot(mapping = aes(x = year, y = value)) +
    geom_point(data = obs, aes(shape = survey)) +
    # geom_point(data = obs, size = .5, alpha = .5)+
    geom_line(data = obs, aes(group = survey, colour = survey), size = .2) +
    geom_ribbon(data = cht, aes(ymin = lower, ymax = upper), fill = 'Grey80', alpha = .3) +
    geom_line(data = cht, linetype = 'dotted', size = 1, colour = 'Grey80') +
    geom_ribbon(data = ind, aes(ymin = lower, ymax = upper), fill = 'Blue', alpha = .1) +
    geom_line(data = ind, linetype = 'solid', size = 1, colour = 'Blue') +
    geom_point(data = influential, aes(x=year, y=value), colour = "mediumseagreen", size=4) +
    #
    scale_y_continuous(limits = c(0,1), breaks = 0:4 * 0.25)+
    scale_x_continuous(limits = c(NA, target_yr))+
    scale_shape_manual(values = 1:n_distinct(obs$survey))+
    labs(x = 'year', y = 'completion rate', colour = 'survey', shape = 'survey')+
    guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
    facet_grid(variable ~ sex, labeller = labeller(.default = label_values)) +
    ggtitle(countrycode::countrycode(target, 'iso3c', 'country.name'))
}

#' Plot ABC Results by Category
#'
#' \code{plt_cat} creates a faceted plot displaying the observed and
#' estimated true completion rates for a given level/sex combination.
#' Each facet presents a different country. The choice of countries can
#' be user provided, or random based on the observations provided.
#'
#' @param mdl_proj Projected values data frame.
#' @param mdl_obs Observed values data frame.
#' @param level Level to plot.
#' @param sexes Sex to plot.
#' @param target_yr Final year to plot.
#' @param countries Vector of countries to plot.
#'
#' @family plotting functions
#'
#' @export
plt_cat <- function(mdl_proj, mdl_obs, level = "prim", sexes = "total",
                    target_yr = 2020, countries = sample(unique(mdl_obs$country), 3)) {

  output <- mdl_proj %>%
    filter(country %in% countries & year <= target_yr & variable == level & sex == sexes) %>%
    mutate(country = countrycode::countrycode(country, 'iso3c', 'country.name'))
  obs <- mdl_obs %>%
    filter(country %in% countries & year <= target_yr & variable == level & sex == sexes) %>%
    mutate(country = countrycode::countrycode(country, 'iso3c', 'country.name'))
  cht <- output %>% filter(series == "projected8")
  ind <- output %>% filter(series == "projected3t5")

  ggplot(mapping = aes(x = year, y = value)) +
    geom_point(data = obs, aes(shape = survey)) +
    # geom_point(data = obs, size = .5, alpha = .5)+
    geom_line(data = obs, aes(group = survey, colour = survey), size = .2) +
    geom_ribbon(data = cht, aes(ymin = lower, ymax = upper), fill = 'Grey80', alpha = .3) +
    geom_line(data = cht, linetype = 'dotted', size = 1, colour = 'Grey80') +
    geom_ribbon(data = ind, aes(ymin = lower, ymax = upper), fill = 'Blue', alpha = .1) +
    geom_line(data = ind, linetype = 'solid', size = 1, colour = 'Blue') +
    #
    scale_y_continuous(limits = c(0,1), breaks = 0:4 * 0.25)+
    scale_x_continuous(limits = c(NA, target_yr))+
    scale_shape_manual(values = 1:n_distinct(obs$survey))+
    labs(x = 'year', y = 'completion rate', colour = 'survey', shape = 'survey') +
    guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
    facet_wrap(~ country)
}

#' Plot Observation Replications
#'
#' \code{mk_repl_plot} creates a plot of 2 observed values from each observed
#' country along with the respective replications generated by the ABC model.
#' The observations are plotted in ascending order of value with the range of
#' the replications around them.
#'
#' @param obs Observed values data frame.
#' @param reps Replications data frame.
#'
#' @family plotting functions
#'
#' @export
mk_repl_plot <- function(obs, reps) {
  obs %>%
    mutate(index = row_number()) %>%
    group_by(country) %>%
    sample_n(2) %>%
    ungroup %>%
    select(country, variable, sex, year, value, index, survey, obsage, recondist, cap_adj) %>%
    rename(target = value) %>%
    dplyr::inner_join(reps, by = c('country', 'year', 'index')) %>%
    mutate(value = pmax(pmin(pnorm(value + cap_adj), 1), 0)) %>%
    arrange(target) %>%
    # mutate(index = ordered(index, levels = unique(index))) %>%
    mutate(index = ordered(interaction(country, index),
                           levels = unique(interaction(country, index)))) %>%
    {ggplot(data = ., aes(x = index, y = value, group = index))+
        # geom_violin(draw_quantiles = c(.1, .5, .9))+
        geom_violin(scale = 'width')+
        # geom_boxplot(outlier.shape = NA, colour = 'grey50')+
        geom_point(aes(y = target), shape = 18, colour = 'red', size = 2)+
        ggthemes::theme_tufte()+
        labs(x = '')+theme(axis.line.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.text.x = element_blank())}
}



#' Plot ABC model bias terms
#'
#' \code{mk_bias_plot} presents survey bias values at each iteration of
#' the model as 90% point ranges. The bias terms are grouped according
#' to the survey source and wave in an effort to present differences,
#' if any, between said sources/waves.
#'
#' @param iters Survey bias iterations.
#'
#' @family plotting functions
#'
#' @export
mk_bias_plot <- function(iters) {
  by_survey <- iters %>%
    group_by(survey, round, level, sex) %>%
    tidybayes::point_interval(value, .width = 0.9) %>%
    select(-.width, -.point, -.interval) %>%
    rename(lower = .lower, upper = .upper, type = round) %>%
    ungroup %>%
    mutate(survey = paste0(survey, level, sex)) %>%
    mutate(survey = forcats::fct_reorder(survey, value))

  by_type <- iters %>%
    select(type=round, iteration, survey, variable, value) %>%
    group_by(type, iteration) %>%
    summarise(value = mean(value, na.rm=TRUE)) %>%
    ungroup %>%
    group_by(type) %>%
    tidybayes::point_interval(value, .width = 0.9) %>%
    select(-.width, -.point, -.interval) %>%
    rename(lower = .lower, upper = .upper) %>%
    ungroup %>%
    select(type, value, lower, upper)

  ggplot() +
    geom_pointrange(data = by_survey,
                    aes(x = survey, y = value, ymin = lower, ymax = upper),#, colour = type),
                    size = .2, fatten = .5)+
    geom_hline(yintercept = 0, size = 0.2)+
    theme_minimal()+
    facet_wrap(~ type, scales = 'free_x')+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())+
    geom_hline(data = by_type, aes(yintercept = value), colour = 'blue', alpha = .4)+
    geom_hline(data = by_type, aes(yintercept = lower), colour = 'blue', alpha = .4, linetype = 'dotted')+
    geom_hline(data = by_type, aes(yintercept = upper), colour = 'blue', alpha = .4, linetype = 'dotted')+
    scale_colour_brewer(type = 'qual', palette = 'Set1')+
    labs(y = 'bias', x = '')
}

#' Plot ABC model parameters
#'
#' \code{mk_par_plot} presents boxplots of each bias/distortion parameter category
#' to illustrate the magnitude of each effect. The parameters in question are:
#' \enumerate{
#'   \item Survey bias
#'   \item Late completion
#'   \item Very late completion
#'   \item Age misreporting
#' }
#'
#' @param df Parameter interval data frame.
#'
#' @family plotting functions
#'
#' @export
mk_par_plot <- function(df) {
  df %>%
    filter(parameter %in% c('beta_s', 'late', 'vlate', 'mult5err')) %>%
    mutate(parameter = factor(parameter, levels = c('beta_s', 'late', 'vlate', 'mult5err'),
                              labels = c('survey bias', 'late completion', 'very late completion', 'age misreporting')),
           level = factor(level, levels = c('prim', 'lsec', 'usec'),
                          labels = c("Primary", "Lower Secondary", "Upper Secondary"))
    ) %>%
    mutate(value = case_when(
      parameter == 'very late completion' ~ value,
      parameter == 'late completion' ~ value/2,
      TRUE ~ value
    )) %>%
    ggplot(data = ., aes(x = parameter, y = value, colour = level))+
    geom_boxplot(width = .3)+
    theme_minimal()+
    geom_hline(yintercept = 0, size = .2)+
    # scale_y_continuous(limits = c(-.5, .5))+
    scale_colour_brewer(type = 'qual', palette = 'Set1')+
    coord_flip() +
    theme(legend.position = "bottom")
}

#' Plot ABC Results by Region
#'
#' \code{mk_reg_plot} plots aggregated ABC model output using either a regional,
#' income level, or mixed aggregation scheme. The outputs present both the
#' official completion rate indicator as well as an estimate of ultimate completion.
#' \code{mk_reg_plot} can be used to produce plots that include forecasts reaching
#' to 2030 though it is worth noting that said forecasts use a mixed within-level
#' and across-levels weighting scheme to recognize the growth across levels is
#' linked in the long term.
#'
#' @param df Regional aggregate data frame.
#' @param which_aggs Type of aggregation (regions, income, or regionsXincome).
#' @param endyr Final year to plot.
#'
#' @family plotting functions
#'
#' @export
mk_reg_plot <- function(df, which_aggs = 'regions', endyr = 2020) {
  plotdat <- df %>%
    filter(year >= 2000 & year <= endyr) %>%
    filter(!(is.na(SDG.region) & is.na(income_group))) %>%
    mutate(income_group = factor(income_group, levels=c("High", "Upper middle", "Lower middle", "Low", "World"))) %>%
    # mutate(income_group = factor(income_group, levels=c("High", "Middle", "Low", "World"))) %>%
    filter(aggregates == 'income' | SDG.region != 'Oceania') %>%
    filter(aggregates == which_aggs) %>%
    mutate(series = case_when(
      series == "projected8" ~ "Ultimate CR",
      series == "projected5" ~ "Upper Bracket CR",
      series == "projected3t5" ~ "CR Indicator",
      TRUE ~ ""),
      series = factor(series, levels = c("CR Indicator", "Upper Bracket CR", "Ultimate CR"))) %>%
    filter(series != "Upper Bracket CR")

  p <- ggplot(data = plotdat, aes(x = year, y = CR))+
    geom_line(aes(linetype = series, colour = level), size = .75)+
    scale_color_brewer(palette = 'Dark2')+
    facet_grid(SDG.region ~ income_group)+
    scale_y_continuous(limits = c(0,1))+
    theme_minimal()

  if (endyr > 2020) {
    p <- p + geom_vline(xintercept = 2020)
  }

  switch(which_aggs,
         regions        = p+facet_wrap(~ SDG.region, ncol = 1),
         income         = p+facet_wrap(~ income_group, ncol = 1),
         regionsXincome = p+facet_wrap(SDG.region ~ income_group) +
           theme(legend.position = "bottom"))
}

#' Plot standard error
#'
#' \code{se_plotter} plots the observed sampling standard errors on
#' the observed scale in a histogram.
#'
#' @param df Input data frame.
#'
#' @family plotting functions
#'
#' @export
se_plotter <- function(df) {
  {ggplot(data = df, aes(x = 100 * se, group = country))+
      geom_histogram(bins = 60, aes(fill = country), alpha = .8, position = 'dodge', na.rm=TRUE)+
      # geom_density(aes(colour = country), bw = .1)+
      scale_x_continuous(breaks = seq(0, 12, 2), limits = c(0, 12))+
      labs(x = 'standard error (p.p.)')+
      theme_minimal()}
}

#' Plot ABC model Rhats
#'
#' \code{mk_rhats_plot} plots the worse case scenario Rhats for the ABC model output, grouped
#' by level and sex.
#'
#' @param rhats_all Model output rhats.
#'
#' @family plotting functions
#'
#' @export
mk_rhats_plot <- function(rhats_all) {
  ggplot(data = rhats_all, aes(x = parameter, y = rhat))+
    theme_bw()+
    geom_hline(yintercept = 1.05, colour = 'Orange', size = 0.3)+
    geom_hline(yintercept = 1.1, colour = 'Red', size = 0.3)+
    geom_segment(aes(xend = parameter, x = parameter, y = rhat, yend = 1), size = 0.3)+
    geom_point(size = 1)+
    coord_flip()+
    scale_y_continuous(limits = c(0.99, 1.11), breaks = c(1, 1.05, 1.1, 1.15))+
    facet_grid(level ~ sex, scales="free_y")+
    theme(panel.grid.major.y = element_blank())
}

#' Plot ABC model trace plots
#'
#' \code{mk_trace_plot} plots the trace plots of the worst case scenario parameters
#' as decided by the rhats diagnostic.
#'
#' @param mcmc ABC model mcmc output.
#' @param rhats_all Model output rhats.
#' @param lvl Level to plot.
#' @param sx Sex to plot.
#'
#' @family plotting functions
#'
#' @export
mk_trace_plot <- function(mcmc, rhats_all, lvl = "prim", sx = "total") {
  pars4trace <-
    rhats_all %>%
    filter(level == lvl, sex == sx) %>%
    pull(parameter)

  bayesplot::mcmc_trace(mcmc, pars = pars4trace, facet_args = list(ncol = 3, nrow = 3),
                        np = bayesplot::nuts_params(mcmc))+theme_minimal()
}
