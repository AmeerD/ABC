# abcR 0.2.0

## Runs on current R / Stan

* Prepped survey files are read with **qs2** instead of **qs**, which has been
  archived on CRAN. `agg_survey()` now reads `data/<survey>_prepped.qs2`; an
  old `.qs` file is reported with conversion instructions instead of being
  misread.
* The ABC independent model now ships with the package (`inst/stan/ABC_indep.stan`,
  path via `abc_model_file()`), written with `array[]` declarations so it compiles
  under Stan >= 2.33, where the old `int x[n]` syntax was removed. The statistical
  model is unchanged.
* `mdl_compile()` compiles the packaged model by default and also accepts a file
  path or Stan code (the old `mdl_compile(readLines(file))` call still works).

## Reproducibility fixes

* `mdl_prep()` gains `end_year`. Previously the latent horizon was
  `year(Sys.Date()) + 5`, so the same data gave models of different size depending
  on the day it was run. Without `end_year` the old behaviour is kept, with a message.
* `mdl_core()` / `mdl_run()` now support several chains in one call (previously
  only one init list was passed, so `nchains > 1` failed). Initial values are
  drawn from a stated seed (`init_seed + chain_nr - 1`, new `mdl_draw_inits()`),
  so one-chain-per-call runs such as `scripts/sample_plans.R` still start each
  chain differently. `seed`, `cores`, `adapt_delta`, `max_treedepth` and
  `refresh` are exposed; the defaults are the previous settings (seed 2019,
  adapt_delta 0.9, max_treedepth 12). The inits used are kept in
  `attr(fit, "abc_init")`.

## New workflow functions

* `abc_example_data()` - the public anonymised test data (`inst/extdata`).
* `add_cap_adj()` - the cap adjustment computed in `scripts/ABCtest.R`, which
  `mdl_prep()`/`mdl_process()` expect in a `cap_adj` column. It stops if a group
  has values beyond both caps, since one shift cannot correct both.
* `split_abc_groups()` - input checks, cap adjustment and the nine level/sex groups.
* `check_abc_input()`, `mdl_check_data()` - input and Stan-data checks, including
  the survey-block ordering the Stan model relies on.
* `mdl_diagnose()`, `mdl_expected_constants()` - Rhat, bulk/tail ESS, divergences,
  tree-depth hits, and a check that zero-variance quantities are exactly the ones
  fixed by the model.
* `abc_run_groups()` - prepare, sample, process and diagnose all groups, with
  saved fits and safe resume.

## Bug fixes

* `get_nsvar()`, `test_abc2()`, `mdl_test()`, `test_cov()` called `pivot_wider()`
  / `replace_na()` without importing them and failed at run time.
* `posterior` is now declared in Imports; unused `readr` removed.
  `R CMD check` passes with no errors, warnings or notes.
