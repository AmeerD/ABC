library(dplyr)
library(abcR)
library(rstan)

load("data/ABCtest.rda")

cap = 0.98
mod = "models/ABC_indep.stan"

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

abc_mod = mdl_compile(readLines(mod))

pt_chain = mdl_run(pt_input, chain_nr = 1, model = abc_mod)
pt_output = mdl_process(pt_input, pt_chain)
pm_chain = mdl_run(pm_input, chain_nr = 1, model = abc_mod)
pm_output = mdl_process(pm_input, pm_chain)
pf_chain = mdl_run(pf_input, chain_nr = 1, model = abc_mod)
pf_output = mdl_process(pf_input, pf_chain)
lt_chain = mdl_run(lt_input, chain_nr = 1, model = abc_mod)
lt_output = mdl_process(lt_input, lt_chain)
lm_chain = mdl_run(lm_input, chain_nr = 1, model = abc_mod)
lm_output = mdl_process(lm_input, lm_chain)
lf_chain = mdl_run(lf_input, chain_nr = 1, model = abc_mod)
lf_output = mdl_process(lf_input, lf_chain)
ut_chain = mdl_run(ut_input, chain_nr = 1, model = abc_mod)
ut_output = mdl_process(ut_input, ut_chain)
um_chain = mdl_run(um_input, chain_nr = 1, model = abc_mod)
um_output = mdl_process(um_input, um_chain)
uf_chain = mdl_run(uf_input, chain_nr = 1, model = abc_mod)
uf_output = mdl_process(uf_input, uf_chain)

output = bind_rows(pt_output, pm_output, pf_output,
                   lt_output, lm_output, lf_output,
                   ut_output, um_output, uf_output)

