library(dplyr)
library(abcR)
library(rstan)

load("data/ABCtest.rda")

cap = 0.98
mod = "models/ABC_indep.stan"

countries = test.dat %>% select(country) %>% distinct %>% pull

caps = test.dat %>%
  group_by(country, variable, sex) %>% 
  summarise(upper = max(value), lower = min(value)) %>%
  ungroup %>%
  mutate(cap_adj = case_when(upper > cap ~ upper - cap,
                             lower < (1-cap) ~ lower - (1-cap),
                             TRUE ~ 0)) %>%
  select(-upper, -lower) %>% 
  identity

pt_input = test.dat %>% filter(variable == "prim" & sex == "total") %>%
  left_join(caps, by = c("country", "variable", "sex")) %>% identity
pm_input = test.dat %>% filter(variable == "prim" & sex == "male") %>%
  left_join(caps, by = c("country", "variable", "sex")) %>% identity
pf_input = test.dat %>% filter(variable == "prim" & sex == "female") %>%
  left_join(caps, by = c("country", "variable", "sex")) %>% identity
lt_input = test.dat %>% filter(variable == "lsec" & sex == "total") %>%
  left_join(caps, by = c("country", "variable", "sex")) %>% identity
lm_input = test.dat %>% filter(variable == "lsec" & sex == "male") %>%
  left_join(caps, by = c("country", "variable", "sex")) %>% identity
lf_input = test.dat %>% filter(variable == "lsec" & sex == "female") %>%
  left_join(caps, by = c("country", "variable", "sex")) %>% identity
ut_input = test.dat %>% filter(variable == "usec" & sex == "total") %>%
  left_join(caps, by = c("country", "variable", "sex")) %>% identity
um_input = test.dat %>% filter(variable == "usec" & sex == "male") %>%
  left_join(caps, by = c("country", "variable", "sex")) %>% identity
uf_input = test.dat %>% filter(variable == "usec" & sex == "female") %>%
  left_join(caps, by = c("country", "variable", "sex")) %>% identity

observed = bind_rows(pt_input, pm_input, pf_input,
                     lt_input, lm_input, lf_input,
                     ut_input, um_input, uf_input)

abc_mod = mdl_compile(readLines(mod))

pt_chain = mdl_run(pt_input, chain_nr = 1, model = abc_mod, nchains=1, nburn=5000, niter=5000, nthin=5)
pt_output = mdl_process(pt_input, pt_chain)
pm_chain = mdl_run(pm_input, chain_nr = 1, model = abc_mod, nchains=1, nburn=5000, niter=5000, nthin=5)
pm_output = mdl_process(pm_input, pm_chain)
pf_chain = mdl_run(pf_input, chain_nr = 1, model = abc_mod, nchains=1, nburn=5000, niter=5000, nthin=5)
pf_output = mdl_process(pf_input, pf_chain)
lt_chain = mdl_run(lt_input, chain_nr = 1, model = abc_mod, nchains=1, nburn=5000, niter=5000, nthin=5)
lt_output = mdl_process(lt_input, lt_chain)
lm_chain = mdl_run(lm_input, chain_nr = 1, model = abc_mod, nchains=1, nburn=5000, niter=5000, nthin=5)
lm_output = mdl_process(lm_input, lm_chain)
lf_chain = mdl_run(lf_input, chain_nr = 1, model = abc_mod, nchains=1, nburn=5000, niter=5000, nthin=5)
lf_output = mdl_process(lf_input, lf_chain)
ut_chain = mdl_run(ut_input, chain_nr = 1, model = abc_mod, nchains=1, nburn=5000, niter=5000, nthin=5)
ut_output = mdl_process(ut_input, ut_chain)
um_chain = mdl_run(um_input, chain_nr = 1, model = abc_mod, nchains=1, nburn=5000, niter=5000, nthin=5)
um_output = mdl_process(um_input, um_chain)
uf_chain = mdl_run(uf_input, chain_nr = 1, model = abc_mod, nchains=1, nburn=5000, niter=5000, nthin=5)
uf_output = mdl_process(uf_input, uf_chain)

output = bind_rows(pt_output, pm_output, pf_output,
                   lt_output, lm_output, lf_output,
                   ut_output, um_output, uf_output)

#A variant of plt_country is needed as the countries in the test dataset have been scrambled and thus 
#do not correspond to actual country codes
label_values <- c(
  "prim" = "Primary",
  "lsec" = "Lower Secondary",
  "usec" = "Upper Secondary",
  "female" = "Female",
  "male" = "Male",
  "total" = "Total"
)

plt_country <- function(mdl_proj, mdl_obs, mdl_influential, target, use_influential = TRUE, target_yr = 2020,
                        levels = c("prim", "lsec", "usec"), sexes = c("female", "male", "total")) {
  output <- mdl_proj %>%
    filter(country == target & year <= target_yr & variable %in% levels & sex %in% sexes) %>%
    #mutate(country = countrycode::countrycode(country, 'iso3c', 'country.name')) %>%
    mutate(variable = factor(variable, levels = levels))
  obs <- mdl_obs %>%
    filter(country == target & year <= target_yr & variable %in% levels & sex %in% sexes) %>%
    #mutate(country = countrycode::countrycode(country, 'iso3c', 'country.name')) %>%
    mutate(variable = factor(variable, levels = levels))
  influential <- mdl_influential %>%
    filter(country == target & year <= target_yr & variable %in% levels & sex %in% sexes) %>%
    #mutate(country = countrycode::countrycode(country, 'iso3c', 'country.name')) %>%
    mutate(variable = factor(variable, levels = levels))
  cht <- output %>% filter(series == "projected8")
  ind <- output %>% filter(series == "projected3t5")
  
  if (!use_influential) {
    influential <- influential %>% anti_join(influential)
  }
  
  ggplot(mapping = aes(x = year, y = value)) +
    geom_point(data = obs, aes(shape = survey)) +
    # geom_point(data = obs, size = .5, alpha = .5)+
    geom_line(data = obs, aes(group = survey, colour = survey), size = .2) +
    geom_ribbon(data = cht, aes(ymin = lower, ymax = upper), fill = 'Grey80', alpha = .3) +
    geom_line(data = cht, linetype = 'dotted', size = 1, colour = 'Grey80') +
    geom_ribbon(data = ind, aes(ymin = lower, ymax = upper), fill = 'Blue', alpha = .1) +
    geom_line(data = ind, linetype = 'solid', size = 1, colour = 'Blue') +
    geom_point(data = influential, aes(x=year, y=value), colour = "mediumseagreen", size=4) +
    #
    scale_y_continuous(limits = c(0,1), breaks = 0:4 * 0.25)+
    scale_x_continuous(limits = c(NA, target_yr))+
    scale_shape_manual(values = 1:n_distinct(obs$survey))+
    labs(x = 'year', y = 'completion rate', colour = 'survey', shape = 'survey')+
    guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
    facet_grid(variable ~ sex, labeller = labeller(.default = label_values))
}

for (c in countries) {
  plt_country(output, observed, observed, c, use_influential = FALSE, target_yr = 2021)
}