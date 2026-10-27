test_that("packaged Stan model parses with current stanc", {
  expect_true(file.exists(abc_model_file()))
  res <- rstan::stanc(file = abc_model_file())
  expect_true(res$status)
})

test_that("example data splits into nine capped groups", {
  inputs <- split_abc_groups(abc_example_data())
  expect_named(inputs, c("pt", "pm", "pf", "lt", "lm", "lf", "ut", "um", "uf"))
  expect_true(all(vapply(inputs, function(d) "cap_adj" %in% names(d), logical(1))))
  # after the shift no value exceeds the cap
  all_vals <- do.call(rbind, inputs)
  expect_true(all(all_vals$value - all_vals$cap_adj <= 0.98 + 1e-12))
})

test_that("mdl_prep uses the explicit horizon and passes structural checks", {
  inputs <- split_abc_groups(abc_example_data())
  dat <- mdl_prep(inputs$lt, end_year = 2031)
  expect_identical(dat$nyears, as.integer(2031 - (min(inputs$lt$year) - 1)))
  expect_silent(mdl_check_data(dat))
  expect_error(mdl_prep(inputs$lt, end_year = 1990), "earlier than")
  expect_message(mdl_prep(inputs$lt), "end_year")
})

test_that("mdl_prep refuses input without cap_adj", {
  df <- abc_example_data()
  df <- df[df$variable == "prim" & df$sex == "total", ]
  expect_error(mdl_prep(df, end_year = 2031), "cap_adj")
})

test_that("inits are reproducible and one per chain", {
  dat <- mdl_prep(split_abc_groups(abc_example_data())$lt, end_year = 2031)
  a <- mdl_draw_inits(dat, nchains = 4, seed = 1)
  b <- mdl_draw_inits(dat, nchains = 4, seed = 1)
  expect_length(a, 4)
  expect_identical(a, b)
})

test_that("input checks reject impossible values", {
  df <- abc_example_data()
  df$value[1] <- 1
  expect_error(check_abc_input(df), "strictly inside")
})

test_that("seeded inits are reproducible and leave the global RNG alone", {
  dat <- list(n_country = 3, nyears = 10)
  set.seed(1); before <- .Random.seed
  a <- mdl_draw_inits(dat, nchains = 2, seed = 5)
  expect_identical(.Random.seed, before)
  expect_identical(a, mdl_draw_inits(dat, nchains = 2, seed = 5))
  expect_false(identical(a[[1]], a[[2]]))
  expect_false(identical(a[1], mdl_draw_inits(dat, nchains = 1, seed = 6)))
})
