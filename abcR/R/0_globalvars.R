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
    "projected3t5", "projected5", "projected8", "recondist3", "refage", "haslate", "obsage", "start",
    "ltype", ".lower", ".upper", ".draw", ".variable", "muzero", "mult5err", "ll", "sigma_s",
    "survey_series", "yhat", "age", "age_adjustment", "cap", "mdk", "mdk_lift", "mdk_mix", "latest",
    "nsurvey", "mad", "mmad", "halfobsage", "weight", "type", "category", "target", "cluster", ".",
    "Region", "SubRegion", "LDC", "LLDC", "SIDS", "Continent", "AU", "GPE", "CAC", "PCFC",
    "Income", "female", "fprop", "male", "mprop", "avg", "mn", "md", "vc", "mu6.ct", "mu7.ct",
    "a3", "a4", "a5", "a6", "a7", "a8", "decile", "sn", "50%", "qval", "lower_0.8", "upper_0.8",
    "lower_0.9", "upper_0.9", "lower_0.95", "upper_0.95", "in80", "in90", "in95", "coverage",
    "mean.chain.1")
  )
