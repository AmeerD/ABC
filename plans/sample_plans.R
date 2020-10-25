
model_plan <- drake_plan(
  trace = TRUE,
  
  rounds = distinct(agg_age, country, year, survey, round),
  
  # subset countries
  country_sample = rounds %>% select(country) %>% distinct() %>% pull()
  ,
  
  # prepare dataset for modelling
  df_4mdl = mk_dataset_4modeling(rounds) %>%
   filter(country %in% country_sample) %>% 
   ungroup() %>%
   mutate(minyear = min(year)) %>%
   group_by(country, variable, sex) %>% 
   mutate(c.count = n(), c.minyear = min(year), backcast = c.minyear - minyear) %>%
   ungroup %>%
   filter(c.count > 3) %>%
   select(-c.count, -minyear, -c.minyear, -backcast) %>%
   group_by(country, survey, variable, sex) %>%
   mutate(s.count = n()) %>%
   filter(s.count != 1) %>%
   ungroup %>%
   select(-s.count) %>% 
   mutate(index = 1:n()) %>%
   identity,
  
  # Construct a table of extreme value adjustments
  cap_lookup = df_4mdl %>%
      group_by(country, variable, sex) %>%
      summarise(upper = max(value), lower = min(value)) %>%
      ungroup %>%
      mutate(cap_adj = case_when(upper > cap ~ upper - cap,
                                 lower < (1-cap) ~ lower - (1-cap),
                                 TRUE ~ 0)) %>%
      select(-upper, -lower) %>%
      identity,
  
  # subset level-gender
  cr_input = target(
   df_4mdl %>% 
     filter(variable == level & sex == sexes) %>% 
     remove_latest(type = dataset) %>%
     left_join(cap_lookup, by=c("country", "variable", "sex")) %>%
     # nest(-country) %>% 
     # sample_n(c_sample_n) %>% 
     # unnest %>% 
     identity
   ,
   transform = cross(level = !! levels2include, sexes = !! sex2include, dataset = !! datatype,
                     .id = c(level, sexes, dataset)) 
  ),
  
  cr_observed = target(
    dplyr::bind_rows(cr_input) %>% distinct %>% select(country, year, survey, variable, sex, value),
    transform = combine(cr_input)
  ),
  
  cr_mod = target(
    #readLines(models[mod]), #
    mdl_compile(readLines(model)),
    hpc = FALSE
  ),
  
  # run model on individual chains
  cr_mcmc_chain = target(
   mdl_run(cr_input, chain_nr = chain, model = cr_mod),
   transform = cross(cr_input, chain = !! seq_len(N_parchain),
                     .id = c(level, sexes, dataset, chain))
  ),
  
  # combine chains for raw multi-chain output
  cr_mcmc = target(
    sflist2stanfit(list(cr_mcmc_chain)),
    transform = combine(cr_mcmc_chain, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset)),
    hpc = FALSE
  ),

  runtime = target(
    data.frame(get_elapsed_time(cr_mcmc)) %>%
      summarise_all(sum) %>%
      mutate(Level = level,
             Sex = sexes,
             Data = dataset,
             Total.Iter = (N_burnin + N_iter) * N_chains * N_parchain,
             Divergences = get_num_divergent(cr_mcmc),
             Exceed.TD = get_num_max_treedepth(cr_mcmc),
             Min.BFMI = min(get_bfmi(cr_mcmc)),
             Min.NEff = min(summary(cr_mcmc)$summary[,"n_eff"], na.rm=TRUE),
             Max.Rhat = max(summary(cr_mcmc)$summary[,"Rhat"], na.rm=TRUE),
             Avg.TD = mean(get_sampler_params(cr_mcmc, inc_warmup=FALSE)[[1]][,"treedepth__"]),
             Avg.LF = mean(get_num_leapfrog_per_iteration(cr_mcmc))) %>%
      select(Level, Sex, everything()) %>%
      mutate(minutes = round((warmup + sample)/60,0)),
    transform = cross(cr_mcmc, .id = c(level, sexes, dataset))
  ),

  rt_summary = target(
    dplyr::bind_rows(runtime) %>% 
      arrange(Level),
    transform = combine(runtime)
  )
)

analysis_plan <- drake_plan(
  trace = TRUE,
  
  rounds = distinct(agg_age, country, year, survey, round),
  
  # subset countries
  country_sample = rounds %>% select(country) %>% distinct() %>% pull()
  ,
  
  # prepare dataset for modelling
  df_4mdl = mk_dataset_4modeling(rounds) %>%
    filter(country %in% country_sample) %>% 
    ungroup() %>%
    mutate(minyear = min(year)) %>%
    group_by(country, variable, sex) %>% 
    mutate(c.count = n(), c.minyear = min(year), backcast = c.minyear - minyear) %>%
    ungroup %>%
    filter(c.count > 3) %>%
    select(-c.count, -minyear, -c.minyear, -backcast) %>%
    group_by(country, survey, variable, sex) %>%
    mutate(s.count = n()) %>%
    filter(s.count != 1) %>%
    ungroup %>%
    select(-s.count) %>% 
    mutate(index = 1:n()) %>%
    identity,
  
  # Construct a table of extreme value adjustments
  cap_lookup = df_4mdl %>%
    group_by(country, variable, sex) %>%
    summarise(upper = max(value), lower = min(value)) %>%
    ungroup %>%
    mutate(cap_adj = case_when(upper > cap ~ upper - cap,
                               lower < (1-cap) ~ lower - (1-cap),
                               TRUE ~ 0)) %>%
    select(-upper, -lower) %>%
    identity,
  
  # subset level-gender
  cr_input = target(
    df_4mdl %>% 
      filter(variable == level & sex == sexes) %>% 
      remove_latest(type = dataset) %>%
      left_join(cap_lookup, by=c("country", "variable", "sex")) %>%
      # nest(-country) %>% 
      # sample_n(c_sample_n) %>% 
      # unnest %>% 
      identity
    ,
    transform = cross(level = !! levels2include, sexes = !! sex2include, dataset = !! datatype,
                      .id = c(level, sexes, dataset)) 
  ),
  
  cr_observed = target(
    dplyr::bind_rows(cr_input) %>% distinct %>% select(country, year, survey, variable, sex, value),
    transform = combine(cr_input)
  ),
  
  cr_mod = target(
    #readLines(models[mod]), #
    mdl_compile(readLines(model)),
    hpc = FALSE
  ),
  
  # run model on individual chains
  cr_mcmc_chain = target(
    mdl_run(cr_input, chain_nr = chain, model = cr_mod),
    transform = cross(cr_input, chain = !! seq_len(N_parchain),
                      .id = c(level, sexes, dataset, chain))
  ),
  
  # combine chains for raw multi-chain output
  cr_mcmc = target(
    sflist2stanfit(list(cr_mcmc_chain)),
    transform = combine(cr_mcmc_chain, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset)),
    hpc = FALSE
  ),
  
  runtime = target(
    data.frame(get_elapsed_time(cr_mcmc)) %>%
      summarise_all(sum) %>%
      mutate(Level = level,
             Sex = sexes,
             Data = dataset,
             Total.Iter = (N_burnin + N_iter) * N_chains * N_parchain,
             Divergences = get_num_divergent(cr_mcmc),
             Exceed.TD = get_num_max_treedepth(cr_mcmc),
             Min.BFMI = min(get_bfmi(cr_mcmc)),
             Min.NEff = min(summary(cr_mcmc)$summary[,"n_eff"], na.rm=TRUE),
             Max.Rhat = max(summary(cr_mcmc)$summary[,"Rhat"], na.rm=TRUE),
             Avg.TD = mean(get_sampler_params(cr_mcmc, inc_warmup=FALSE)[[1]][,"treedepth__"]),
             Avg.LF = mean(get_num_leapfrog_per_iteration(cr_mcmc))) %>%
      select(Level, Sex, everything()) %>%
      mutate(minutes = round((warmup + sample)/60,0)),
    transform = cross(cr_mcmc, .id = c(level, sexes, dataset))
  ),
  
  rt_summary = target(
    dplyr::bind_rows(runtime) %>% 
      arrange(Level),
    transform = combine(runtime)
  ),
  
  # Process model output
  cr_output = target(
   mdl_process(cr_input, cr_mcmc),
   transform = combine(cr_input, cr_mcmc, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset))
  ),

  # # Combine model output
  cr = target(
   dplyr::bind_rows(cr_output),
   transform = combine(cr_output, .by = c(dataset))
  ),
  
  cr_repl = target(
    get_yrepl(cr_input, cr_mcmc),
    transform = combine(cr_input, cr_mcmc, .by=c(level, sexes, dataset), .id=c(level, sexes, dataset))
  ),

  cr_drift = target(
   get_drift(cr_input, cr_mcmc),
   transform = combine(cr_input, cr_mcmc, .id = c(level, sexes, dataset), .by=c(level, sexes, dataset))
  ),

  cr_drift2 = target(
    cr_drift %>% mutate(Level = level, Sex = sexes),
    transform = cross(cr_drift, .id = c(level, sexes, dataset))
  ),

  cr_drift_all = target(
    dplyr::bind_rows(cr_drift2),
    transform = combine(cr_drift2, .by = c(dataset))
  ),
  
  cr_rhat = target(
    get_rhats(cr_mcmc, level, sexes),
    transform = map(cr_mcmc, .id = c(level, sexes, dataset))
  ),
  
  cr_rhat_all = target(
    dplyr::bind_rows(cr_rhat),
    transform = combine(cr_rhat, .by = c(dataset))
  ),
  
  cr_pars = target(
    get_pars(cr_input, cr_mcmc),
    transform = combine(cr_input, cr_mcmc, .id = c(level, sexes, dataset), .by = c(level, sexes, dataset))
  ),
  
  cr_bias = target(
    get_bias(cr_input, cr_mcmc, rounds),
    transform = combine(cr_input, cr_mcmc, .id = c(level, sexes, dataset), .by = c(level, sexes, dataset))
  ),
  
  cr_bias_all = target(
    dplyr::bind_rows(cr_bias),
    transform = combine(cr_bias, .by = c(dataset))
  ),
  
  cr_pars2 = target(
    dplyr::bind_rows(cr_pars),
    transform = combine(cr_pars, .by = c(sexes, dataset))
  ),
  
  cr_loo = target(
   loo_output(cr_mcmc),
   transform = combine(cr_mcmc, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset)),
   hpc = FALSE
  ),

  cr_loo_influential = target(
   loo_influential(cr_input, cr_loo),
   transform = combine(cr_input, cr_loo, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset)),
   hpc = FALSE
  ),
  
  cr_influential = target(
   dplyr::bind_rows(cr_loo_influential),
   transform = combine(cr_loo_influential, .by = c(dataset))
  ),
  
  #TESTING PHASE
  
  cr_testset = target(
    df_4mdl %>%
      filter(variable == level & sex == sexes) %>%
      dplyr::anti_join(cr_input, by = "index"),
    transform = map(cr_input, .id = c(level, sexes, dataset))
  ),

  cr_test = target(
    test_abc(cr_input, cr_testset, cr_mcmc),
    transform = combine(cr_input, cr_testset, cr_mcmc, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset)),
    hpc = FALSE
  ),

  cr_test_all = target(
    dplyr::bind_rows(cr_test),
    transform = combine(cr_test, .by = c(dataset), .id = c(dataset))
  ),
  
  alt1 = target(
    alt1_smpl2(cr_input, cr_testset),
    transform = combine(cr_input, cr_testset, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset))
  ),

  alt1_all = target(
    dplyr::bind_rows(alt1),
    transform = combine(alt1, .by = c(dataset), .id = c(dataset))
  ),

  alt2 = target(
    alt2_flat2(cr_input, cr_testset),
    transform = combine(cr_input, cr_testset, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset))
  ),

  alt2_all = target(
    dplyr::bind_rows(alt2),
    transform = combine(alt2, .by = c(dataset), .id = c(dataset))
  ),
  
  cr_testcr = target(
    test_abc2(cr_input, cr_testset, cr_mcmc),
    transform = combine(cr_input, cr_testset, cr_mcmc, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset)),
    hpc = FALSE
  ),
  
  cr_testcr_all = target(
    dplyr::bind_rows(cr_testcr),
    transform = combine(cr_testcr, .by = c(dataset), .id = c(dataset))
  ),
  
  alt1cr = target(
    alt1_smpl3(cr_input, cr_testset),
    transform = combine(cr_input, cr_testset, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset))
  ),
  
  alt1cr_all = target(
    dplyr::bind_rows(alt1cr),
    transform = combine(alt1cr, .by = c(dataset), .id = c(dataset))
  ),
  
  alt2cr = target(
    alt2_flat3(cr_input, cr_testset),
    transform = combine(cr_input, cr_testset, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset))
  ),
  
  alt2cr_all = target(
    dplyr::bind_rows(alt2cr),
    transform = combine(alt2cr, .by = c(dataset), .id = c(dataset))
  ),
  
  alt3cr = target(
    alt3_ltst3(cr_input, cr_testset),
    transform = combine(cr_input, cr_testset, .by = c(level, sexes, dataset), .id = c(level, sexes, dataset))
  ),
  
  alt3cr_all = target(
    dplyr::bind_rows(alt3cr),
    transform = combine(alt3cr, .by = c(dataset), .id = c(dataset))
  ),
  
  test_summary = target(
    dplyr::bind_rows(cr_test_all, alt1_all, alt2_all),
    transform = combine(cr_test_all, alt1_all, alt2_all, .by = c(dataset))
  ),
  
  crtest_summary = target(
    dplyr::bind_rows(cr_testcr_all, alt1cr_all, alt2cr_all, alt3cr_all),
    transform = combine(cr_testcr_all, alt1cr_all, alt2cr_all, alt3cr_all, .by=c(dataset))
  )
)

report_plan <- drake_plan(
  #Input plots
  retro_plot = mk_input_plot(cr_input_prim_total_full, c("BDI","SSD")),
  se_plot = se_plotter(df_4mdl %>% filter(country %in% c("BLZ", "MLI"))),
  mult5_plot = mk_input_survey_plot(cr_input_prim_total_full, c("NGADHS2003", "NGAMICS2016")),
  
  #Output plots
  late_plot = plt_cat(cr_full, cr_observed, level="prim", sexes="total", countries=c("MWI")) +
    theme_minimal(),
  prim_country_results = plt_cat(cr_full, cr_observed, level = "prim", sexes = "total",
                                 countries = c("BDI", "HND", "RWA")) +
    theme_minimal(base_size=9) +
    theme(legend.position = "none"),
  lsec_country_results = plt_cat(cr_full, cr_observed, level = "lsec", sexes = "total",
                                 countries = c("BTN", "DNK", "KEN")) +
    theme_minimal(base_size=9) +
    theme(legend.position = "none"),
  usec_country_results = plt_cat(cr_full, cr_observed, level = "usec", sexes = "total",
                                 countries = c("ARM", "KGZ", "LCA")) +
    theme_minimal(base_size=9) +
    theme(legend.position = "none"),
  
  #Regional Aggregation
  cr_reg = combine_regs(cr_full),
  cr_reg2 = target(
    cr_reg %>% filter(sex == s) %>% select(-sex),
    transform = map(s = !!sex2include, .id=c(s))
  ),
  cr_proj = target(
    proj2030(cr_reg2),
    transform = map(cr_reg2, .id=c(s))
  ),
  reg_plot = target(
    mk_reg_plot(cr_reg2, aggregates),
    transform = cross(cr_reg2, aggregates = !! c('regions', 'income', 'regionsXincome'),
                      .id = c(s, aggregates))
  ),
  proj_plot = target(
    mk_reg_plot(cr_proj, aggregates, endyr = 2030),
    transform = cross(cr_proj, aggregates = !! c('regions', 'income', 'regionsXincome'),
                      .id = c(s, aggregates))
  ),
  
  #Gamma Dispersion Plot
  gamma_countries = cr_drift2_prim_total_full %>% 
    arrange(value) %>%
    mutate(quintile = ntile(value, 5)) %>%
    group_by(quintile) %>%
    sample_n(1) %>%
    ungroup %>%
    select(country, quintile),
  
  gamma_data = cr_mcmc_prim_total_full %>%
    tidybayes::recover_types(cr_input_prim_total_full) %>%
    tidybayes::gather_draws(drift[country]) %>%
    ungroup %>%
    select(-.chain, -.iteration, -.variable) %>%
    rename(iteration = .draw, value = .value) %>%
    inner_join(gamma_countries, by = "country"),
  
  gamma_plot = gamma_data %>%
    mutate(country = countrycode::countrycode(country, 'iso3c', 'country.name')) %>%
    ggplot() +
    geom_density(aes(x=value, group=country, colour=country)) +
    stat_function(data=data.frame(x=seq(-0.25, 0.75, by=0.01)), aes(x), fun=emg::demg, 
                  args=c(mu=0, sigma=0.025, lambda=2), colour="black", n=1001) +
    scale_y_continuous(trans='sqrt') +
    theme_minimal() +
    theme(legend.position = "bottom"),
  
  #Diagnostics Plots
  repl_plot_prim = mk_repl_plot(cr_input_prim_total_full, cr_repl_prim_total_full),
  repl_plot_lsec = mk_repl_plot(cr_input_lsec_total_full, cr_repl_lsec_total_full),
  repl_plot_usec = mk_repl_plot(cr_input_usec_total_full, cr_repl_usec_total_full),
  pars_plot = mk_par_plot(cr_pars2_total_full),
  bias_plot = mk_bias_plot(cr_bias_all_full),
  rhats_plot = mk_rhats_plot(cr_rhat_all_full %>%
                              group_by(level, sex) %>%
                              arrange(rhat) %>%
                              mutate(parameter = row_number())) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
  
  trace_plot = mk_trace_plot(cr_mcmc_prim_total_full, cr_rhat_all_full),
  
  countries_sorted = sort(country_sample),

  results = write.csv(cr_full, file_out("./output/cr_full.csv")),
  
  plt = target(
   plt_country(cr_full, cr_observed, cr_influential_full, countries_sorted),
   dynamic = map(sort(countries_sorted))
  ),

  report = {
    rmarkdown::render(
      knitr_in('ABCplots.Rmd'),
      output_file = file_out('ABCplots.html'),
      quiet = TRUE
    )
    #file_out('output/ABCplots.html')
    },
  
  paper = rmarkdown::render(
    knitr_in('./output/CR technical paper v2.Rmd'),
    output_file = file_out('./output/CR-technical-paper-v2.pdf'),
    quiet = TRUE
  )
)

output_plan <- bind_plans(analysis_plan, report_plan)