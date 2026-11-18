# Tests for add_cap_adj() and split_abc_groups(). No Stan sampling is needed.
local({
  suppressMessages(library(dplyr))
  library(abcR)

  df <- data.frame(
    country  = rep(c("AAA", "BBB", "CCC"), each = 4),
    variable = "prim",
    sex      = rep(c("total", "male"), times = 6),
    value    = c(0.990, 0.50, 0.970, 0.60,   # AAA total: max 0.99 -> shift down
                 0.010, 0.40, 0.300, 0.70,   # BBB total: min 0.01 -> shift up
                 0.500, 0.30, 0.600, 0.40),  # CCC: inside the caps
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(add_cap_adj(df))

  # rows and existing columns are untouched, one column is added
  stopifnot(identical(out[names(df)], df), identical(setdiff(names(out), names(df)), "cap_adj"))
  adj <- function(country, sex) unique(out$cap_adj[out$country == country & out$sex == sex])
  stopifnot(isTRUE(all.equal(adj("AAA", "total"), 0.99 - 0.98)))
  stopifnot(isTRUE(all.equal(adj("BBB", "total"), 0.01 - 0.02)))
  stopifnot(identical(adj("CCC", "total"), 0), identical(adj("AAA", "male"), 0))
  # after the shift every value lies inside the caps
  shifted <- out$value - out$cap_adj
  stopifnot(all(shifted <= 0.98 + 1e-12), all(shifted >= 0.02 - 1e-12))

  # an existing cap_adj column is replaced, not duplicated
  again <- suppressMessages(add_cap_adj(out))
  stopifnot(identical(again, out))

  # values on both sides of the caps: shifted down, with a warning
  both <- data.frame(country = "DDD", variable = "prim", sex = "total",
                     value = c(0.995, 0.005), stringsAsFactors = FALSE)
  warned <- FALSE
  res <- withCallingHandlers(
    suppressMessages(add_cap_adj(both)),
    warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
  stopifnot(warned, isTRUE(all.equal(unique(res$cap_adj), 0.995 - 0.98)))

  # groups: named by level and sex, empty groups left out, rows preserved
  groups <- suppressMessages(split_abc_groups(df))
  stopifnot(identical(names(groups), c("pt", "pm")))
  expected <- out[out$sex == "total", ]
  rownames(expected) <- NULL
  got <- as.data.frame(groups$pt)
  rownames(got) <- NULL
  stopifnot(identical(got, expected))
  stopifnot(sum(vapply(groups, nrow, integer(1))) == nrow(df))
})
