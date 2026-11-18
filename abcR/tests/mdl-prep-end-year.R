# Tests for the `end_year` argument of mdl_prep(). No Stan sampling is needed.
local({
  suppressMessages(library(dplyr))
  library(abcR)

  # a small synthetic input for one level/sex group, in the layout of
  # data/ABCtest.rda: two countries, two surveys each, eight cohorts per survey
  df <- expand.grid(k = 0:7, obsyear = c(2005, 2012), country = c("AAA", "BBB"),
                    stringsAsFactors = FALSE) %>%
    mutate(survey = paste0(country, obsyear), year = obsyear - k,
           variable = "prim", sex = "total",
           value = 0.75 - 0.02 * k + ifelse(country == "AAA", 0.1, 0),
           se = 0.02, obsage = 2 - k, recondist = pmax(0, k - 2),
           truage5mlt = 0L, cap_adj = 0) %>%
    select(-k)
  baseyear <- min(df$year) - 1

  # default: current year plus five, exactly as before
  default <- mdl_prep(df)
  stopifnot(identical(default$nyears,
                      lubridate::year(Sys.Date()) + 5 - baseyear))

  # an explicit year fixes the horizon and changes nothing else
  fixed <- mdl_prep(df, end_year = 2030)
  stopifnot(identical(fixed$nyears, 2030 - baseyear))
  stopifnot(identical(names(fixed), names(default)))
  stopifnot(identical(fixed[names(fixed) != "nyears"],
                      default[names(default) != "nyears"]))

  # invalid values are rejected
  fails <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
  stopifnot(fails(mdl_prep(df, end_year = 2005)))
  stopifnot(fails(mdl_prep(df, end_year = 2030.5)))
  stopifnot(fails(mdl_prep(df, end_year = c(2030, 2031))))
})
