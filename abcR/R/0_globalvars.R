#List of dplyr column names used in the package
# These must be added as global variables to address the
# "no visible binding" R CMD check note
utils::globalVariables(
  c("variable", "level", "country", "survey", "series", "sex", "value", "lower", "upper",
    "model", "data", "year", "pred", "iteration", "diff_sq", "recondist", "obsyear",
    "truage5mlt", "regions", "SDG.region", "income_group", "beta_s", ".iteration", ".chain",
    ".varaible", ".value", ".point", ".interval", ".width", "drift", "late", "vlate", "mu_ct",
    "n", "yrepl", "Pareto.k", "modelfile", "se", "se_q", "se_q2", "cap_adj", "lmean", "irr",
    "late_id", "m", "se_temp", "raw_value", "was_capped", "agg_age", "agg_age_gender", "refages",
    "parameter", "aggregates", "CR", "index", "rhat", "mse", "mu3.ct", "mu4.ct", "mu5.ct", "mu8.ct",
    "projected3t5", "projected5", "projected8", "recondist3", "refage", "haslate", "obsage",
    "ltype", ".lower", ".upper", ".draw", ".variable", "muzero", "mult5err", "ll", "sigma_s",
    "survey_series", "yhat", "age", "age_adjustment", "cap", "mdk", "mdk_lift", "mdk_mix", "latest",
    "nsurvey", "mad", "mmad", "halfobsage", "weight", "type")
  )
