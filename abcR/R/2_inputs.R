# MODEL INPUT HELPERS ------------------------------------------------------

#' Add the cap adjustment to ABC input data.
#'
#' Observed completion rates very close to 0 or 1 have unstable probit
#' transforms. For every country, level and sex, if the largest value exceeds
#' \code{cap} all values of that group are shifted down by the excess; otherwise,
#' if the smallest value is below \code{1 - cap}, they are shifted up by the
#' shortfall. The shift is stored in the column \code{cap_adj}, which
#' \code{\link{mdl_prep}} subtracts before the probit transformation and
#' \code{\link{mdl_process}} adds back. This is the calculation carried out at
#' the top of \code{scripts/ABCtest.R}.
#'
#' A group with values on both sides (above \code{cap} and below
#' \code{1 - cap}) cannot be brought inside the caps by a single shift; it is
#' shifted down, as in the script, and a warning names the group.
#'
#' @param df Input data frame with columns \code{country}, \code{variable},
#'   \code{sex} and \code{value}.
#' @param cap Extreme value cap.
#'
#' @return \code{df} with the column \code{cap_adj} added (or replaced).
#'
#' @seealso \code{\link{split_abc_groups}}
#' @export
add_cap_adj <- function(df, cap = 0.98) {
  caps <- df %>%
    group_by(country, variable, sex) %>%
    summarise(upper = max(value), lower = min(value)) %>%
    ungroup
  both <- caps$upper > cap & caps$lower < (1 - cap)
  if (any(both)) {
    warning("Values both above ", cap, " and below ", 1 - cap, " in: ",
            paste(caps$country[both], caps$variable[both], caps$sex[both],
                  collapse = "; "),
            ". These groups are shifted down only.", call. = FALSE)
  }
  caps <- caps %>%
    mutate(cap_adj = case_when(upper > cap ~ upper - cap,
                               lower < (1 - cap) ~ lower - (1 - cap),
                               TRUE ~ 0)) %>%
    select(-upper, -lower)
  if ("cap_adj" %in% names(df)) {
    df <- df %>% select(-cap_adj)
  }
  df %>%
    left_join(caps, by = c("country", "variable", "sex"))
}

#' Split ABC input data into the level and sex groups that are modelled.
#'
#' The ABC model is run separately for each education level and sex.
#' \code{split_abc_groups} adds the cap adjustment with
#' \code{\link{add_cap_adj}} and returns one data frame per group, ready for
#' \code{\link{mdl_run}}. It replaces the nine \code{filter()} and
#' \code{left_join()} calls in \code{scripts/ABCtest.R}.
#'
#' @param df Input data frame with columns \code{country}, \code{variable}
#'   (\code{"prim"}, \code{"lsec"} or \code{"usec"}), \code{sex}
#'   (\code{"total"}, \code{"male"} or \code{"female"}) and \code{value}.
#' @param cap Extreme value cap, passed to \code{\link{add_cap_adj}}.
#'
#' @return A named list of data frames. The names combine the first letter of
#'   the level and of the sex: \code{pt}, \code{pm}, \code{pf}, \code{lt},
#'   \code{lm}, \code{lf}, \code{ut}, \code{um} and \code{uf}. Groups without
#'   data are left out.
#'
#' @export
split_abc_groups <- function(df, cap = 0.98) {
  df <- add_cap_adj(df, cap = cap)
  levels <- c("prim", "lsec", "usec")
  sexes  <- c("total", "male", "female")
  out <- list()
  for (l in levels) {
    for (s in sexes) {
      group <- df %>% filter(variable == l & sex == s)
      if (nrow(group) > 0) {
        out[[paste0(substr(l, 1, 1), substr(s, 1, 1))]] <- group
      }
    }
  }
  out
}
