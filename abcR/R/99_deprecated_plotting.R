#' Core Plotting Function for ABC Results (Deprecated)
#'
#' @param mdl_out Model output (Processed by deprecated mdl_process).
#' @param country_selection Vector of country codes to plot.
#' @param target_yr Final year to plot, defaults to 2020.
#'
#' @export
mdl_plot <- function(mdl_out, country_selection, target_yr = 2020) {
  .Deprecated("plt_country")
  output <- mdl_out$crs
  d <- purrr::map(output, .f = ~ {dplyr::filter(., country %in% country_selection & year <= target_yr) %>%
      mutate(country = countrycode::countrycode(country, 'iso3c', 'country.name'))})
  obs <- d$observed
  cht <- d$projected8
  ind <- d$projected3t5

  ggplot(mapping = aes(x = year, y = value))+
    geom_point(data  = obs, aes(shape = survey))+
    # geom_point(data = obs, size = .5, alpha = .5)+
    geom_line(data   = obs, aes(group = survey, colour = survey), size = .2)+
    geom_ribbon(data = cht, aes(ymin = lower, ymax = upper), fill = 'Grey80', alpha = .3)+
    geom_line(data   = cht, linetype = 'dotted', size = 1, colour = 'Grey80')+
    geom_ribbon(data = ind, aes(ymin = lower, ymax = upper), fill = 'Blue', alpha = .1)+
    geom_line(data   = ind, linetype = 'solid', size = 1, colour = 'Blue')+
    #
    scale_y_continuous(limits = c(0,1), breaks = 0:4 * 0.25)+
    scale_x_continuous(limits = c(NA, target_yr))+
    scale_shape_manual(values = 1:n_distinct(obs$survey))+
    labs(x = 'year', y = 'completion rate', colour = 'survey', shape = 'survey')+
    guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))+
    facet_wrap(~ country, nrow = 1)+
    theme_minimal()
}

#' Plot ABC Results by Gender (Deprecated)
#'
#' @param mdl_out_g Model output (Processed by deprecated mdl_process).
#' @param ... Pass through parameters.
#'
#' @export
mdl_plot_gender <- function(mdl_out_g, ...) {
  mdl_plot(mdl_out_g, ...)+facet_wrap(~ gender, nrow = 1)
}

#' Plot alternate model outputs (Deprecated)
#'
#' @param df Simple model output data frame.
#' @param country_selection Vector of countries to plot.
#'
#' @export
smpl_plot <- function(df, country_selection) {
  ggplot(data = filter(df, country %in% country_selection),
         aes(x = year))+
    geom_point(aes(y = value))+
    geom_line(aes(y = pred))+
    geom_point(aes(y = pred), shape = 5)+
    scale_y_continuous(limits = c(0,1), breaks = 0:4 * .25)+
    labs(x = 'year', y = 'primary completion rate')+
    facet_wrap(level ~ country, ncol = 1)+
    theme_minimal()
}

#' Observation replication plot (Deprecated)
#'
#' @param tgt Observed values data frame.
#' @param df Replications data frame.
#'
#' @export
mk_repl_plot2 <- function(tgt, df) {
  tgt %>%
    filter(recondist %in% c(3, 8, 18)) %>%
    group_by(interaction(obsage, recondist)) %>%
    # distinct(country, index, .keep_all = TRUE) %>%
    # sample_n(., min(n_groups(.), 10)) %>%
    ungroup %>%
    rename(target = value) %>%
    dplyr::inner_join(df, by = c('country', 'year', 'index')) %>%
    mutate(value = pnorm(value)) %>%
    arrange(target) %>%
    # mutate(index = ordered(index, levels = unique(index))) %>%
    mutate(index = ordered(interaction(country, index),
                           levels = unique(interaction(country, index)))) %>%
    {ggplot(data = ., aes(x = index, y = value, group = index))+
        # geom_violin(draw_quantiles = c(.1, .5, .9))+
        geom_violin(scale = 'width')+
        # geom_boxplot(outlier.shape = NA, colour = 'grey50')+
        geom_point(aes(y = target), shape = 18, colour = 'red', size = 3)+
        ggthemes::theme_tufte()+
        labs(x = '')+theme(axis.line.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.text.x = element_blank())}
}

#' Plot input backprojection (Deprecated)
#'
#' @param df Input data frame.
#' @param cs Countries to plot.
#' @param lvl Level to plot.
#' @param sx Sex to plot.
#'
#' @export
bkprj_plot <- function(df, cs, lvl = 'prim', sx = 'total') {
  if (is.null(cs) | !cs %in% unique(df$country)) {
    cs <- unique(df$country)
  }
  plotdat <- filter(df, country %in% cs & variable == lvl & sex == sx)
  blacklisted <- blacklist(plotdat) %>%
    mutate(
      obsage = factor(obsage + 3),
      country = factor(country))
  whitelisted <- whitelist(plotdat) %>%
    mutate(
      obsage = factor(obsage + 3),
      country = factor(country))

  plotdat %>%
    {ggplot(data = ., aes(x = year, y = value, group = survey))+
        geom_line(aes(colour = survey))+
        geom_point(data = filter(plotdat, truage5mlt == 1), shape = 3, size = 2)+
        geom_point(data = whitelisted, aes(alpha = 1 - recondist, shape = obsage), size = 2)+
        scale_y_continuous(limits = c(0,1))+
        facet_wrap(~ country)+
        labs(shape = 'age observed\n (above nominal graduation age)', alpha = 'backprojection distance')+
        guides(
          # colour = guide_legend(ncol = 2, nrow = 3),
          shape = guide_legend(direction = 'horizontal', title.position = 'top', label.position = 'bottom'),
          alpha = guide_legend(direction = 'horizontal', title.position = 'top', label.position = 'bottom'))+
        theme_minimal()+
        geom_point(data = blacklisted, colour = 'Red', shape = 25, size = 2, fill = 'Red', alpha = 0.5)}
}

#' Plot long-term drift (Deprecated)
#'
#' @param df Model output (Processed by deprecated mdl_process).
#'
#' @export
mk_drift_plot <- function(df) {
  drifts <-
    df$pars_iter$drift %>%
    {semi_join(
      .,
      group_by(., country) %>%
        summarise(value = mean(value)) %>%
        mutate(q = cut(value, c(-Inf, quantile(value), Inf), labels = FALSE)) %>%
        group_by(q) %>%
        sample_n(1) %>%
        ungroup,
      by = 'country')}

  drift_prior <- data.frame(x = rnorm(1e7, sd = 1/100, mean = rexp(1e7, rate = 1)))

  ggplot(data = drift_prior, aes(x = x))+
    geom_density(adjust = .2, aes(y = ..scaled.. * 100), n = 8192)+
    theme_minimal()+
    coord_cartesian(xlim = c(-.1, 1.5))+
    scale_x_continuous(breaks = seq(-.1, 1.5, .1))+
    geom_density(data = gs, adjust = 5, n = 8192,
                 aes(group = country, fill = factor(country),
                     x = value, y = ..scaled.. * 100), alpha = 0.5, size = .1)+
    labs(y = 'density (scaled)', x = 'drift', fill = 'country')
}
