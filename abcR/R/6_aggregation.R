#' Regional aggregation
#'
#' Aggregates ABC data based on population weights by country and level.
#'
#' @param df Input data frame.
#' @param pops Population weights.
reg_aggs <- function(df, pops) {
  df %>%
    dplyr::left_join(pops, by = c('country', 'year', 'series', 'sex', 'level')) %>%
    group_by_at(.vars = vars(-value, -country, -weight)) %>%
    dplyr::summarise(CR = round(weighted.mean(value, w = weight, na.rm = TRUE), 3)) %>%
    ungroup
}

#' Complete regional aggregation
#'
#' \code{combine_regs} aggregates completion rate estimates based on SDG region,
#' income group, and SDG region x income group using population weights. The
#' results of each type of aggregation are then binded into a single data frame.
#'
#' @param df Input data frame.
#' @param regions SDG regions table.
#' @param pops Population weights by country, year, series, sex, and level.
#'
#' @return Data frame summarising aggregations.
#' @export
combine_regs <- function(df, regions, pops) {
  poptemp <- pops %>%
    mutate(series = case_when(
      series == "cr" ~ "projected3t5",
      series == "a5" ~ "projected5",
      series == "a8" ~ "projected8",
      TRUE ~ series
    ))

  df %>%
    filter(series != 'observed') %>%
    select(-lower, -upper) %>%
    mutate(level = case_when(
      variable == "prim" ~ "primary",
      variable == "lsec" ~ "lower secondary",
      variable == "usec" ~ "upper secondary",
      variable == "four" ~ "four years",
      TRUE ~ ""
    )) %>%
    select(-variable) %>%
    {bind_rows(
      dplyr::left_join(., regions, by = 'country')
      ,
      .
    )} %>%
    mutate_at(vars(-country, -year, -series, -value, -sex, -level),
              function(x) ifelse(is.na(x), "World", x)) %>%
    {bind_rows(
      select(., -SDG.region, -LDC, -LLDC, -SIDS, -continent, -AU, -GPE, -FCS) %>%
        {bind_rows(
          .,
          filter(., income_group %in% c("Low", "Lower middle")) %>% mutate(income_group = "Low and lower middle"),
          filter(., income_group %in% c("Lower middle", "Upper middle")) %>% mutate(income_group = "Lower and Upper middle")
        )} %>%
        reg_aggs(poptemp) %>%
        mutate(aggregates = 'income'),
      select(., -income_group, -LDC, -LLDC, -SIDS, -continent, -AU, -GPE, -FCS) %>%
        reg_aggs(poptemp) %>%
        mutate(aggregates = 'regions'),
      select(., -LDC, -LLDC, -SIDS, -continent, -AU, -GPE, -FCS) %>%
        reg_aggs(poptemp) %>%
        mutate(aggregates = 'regionsXincome'),
      select(., -SDG.region, -income_group, -LDC, -LLDC, -SIDS, -AU, -GPE, -FCS) %>%
        reg_aggs(poptemp) %>%
        mutate(aggregates = 'continent'),
      select(., -SDG.region, -income_group, -LLDC, -SIDS, -continent, -AU, -GPE, -FCS) %>%
        reg_aggs(poptemp) %>%
        mutate(aggregates = 'LDC'),
      select(., -SDG.region, -income_group, -LDC, -SIDS, -continent, -AU, -GPE, -FCS) %>%
        reg_aggs(poptemp) %>%
        mutate(aggregates = 'LLDC'),
      select(., -SDG.region, -income_group, -LDC, -LLDC, -continent, -AU, -GPE, -FCS) %>%
        reg_aggs(poptemp) %>%
        mutate(aggregates = 'SIDS'),
      select(., -SDG.region, -income_group, -LDC, -LLDC, -SIDS, -continent, -GPE, -FCS) %>%
        reg_aggs(poptemp) %>%
        mutate(aggregates = 'AU'),
      select(., -SDG.region, -income_group, -LDC, -LLDC, -SIDS, -continent, -AU, -FCS) %>%
        reg_aggs(poptemp) %>%
        mutate(aggregates = 'GPE'),
      select(., -SDG.region, -income_group, -LDC, -LLDC, -SIDS, -continent, -GPE, -AU) %>%
        reg_aggs(poptemp) %>%
        mutate(aggregates = 'FCS')
    )} %>%
    order_levels
}

#' Projections to 2030
#'
#' \code{proj2030} projects aggregated completion rate estimates forward to the
#' year 2030. Projections use a weighted average of within-level and across-level
#' long-term drifts with weights shifting linearly from the former to the latter
#' across years. Such a decision is motivated by the understanding that in the
#' longer term, there is likely to interactions between the levels that is not
#' captured in the short term by the independent modeling. For example, large
#' growth in primary education may eventually translate into greater growth in
#' the secondary levels in the long term.
#'
#' @param df Aggregated data frame.
#' @param horizon Projection horizon (defaults to 10).
#'
#' @return An extended data frame formatted for plotting that includes the projected values.
#'
#' @export
proj2030 <- function(df, horizon = 10) {
  growth <- df %>%
    filter(year >= 2010) %>%
    group_by(aggregates, income_group, SDG.region, series, level) %>%
    summarise(mdk = mean(diff(pmax(-3, pmin(3, qnorm(CR)))), na.rm = TRUE)) %>%
    arrange(level) %>%
    # mutate(mdk_lift = map2_dbl(mdk, cummax(mdk), ~ p * .x + (1-p) * .y)) %>%
    mutate(mdk_lift = mean(mdk)) %>%
    arrange(aggregates, income_group, SDG.region, series, level)

  proj <- df %>%
    group_by(aggregates, income_group, SDG.region, series, level) %>%
    filter(year == 2020) %>%
    left_join(growth) %>%
    nest(year, CR, mdk, mdk_lift) %>%
    mutate(data = purrr::map(data, ~ data.frame(
      CR = .$CR, mdk = .$mdk, mdk_lift = .$mdk_lift,
      year = .$year + 0:horizon ,
      weight = pmin(1, {0:horizon}/horizon)
      #weight = 0
    ))) %>%
    mutate(data = purrr::map(data, ~ mutate(.x,
                                            mdk_mix = (1-weight) * mdk + weight * mdk_lift))) %>%
    mutate(data = purrr::map(data, ~ mutate(.x,
                                            CR = pnorm(qnorm(CR) + c(0, cumsum(head(mdk_mix, -1))))))) %>%
    unnest %>%
    select(-mdk, -mdk_lift, -mdk_mix, -weight) %>%
    filter(year > 2020) %>%
    ungroup

  df %>%
    filter(year <= 2020) %>%
    bind_rows(proj)
}
