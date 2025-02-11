
### Bayesian analysis of beach cohort survey data - AGI outcome

pacman::p_load(
  Matrix,
  tidyverse,  
  rstatix,
  janitor,
  brms,
  tidybayes, 
  bayesplot,
  marginaleffects,
  cmdstanr
)

# Start with AGI outcome, first build model with no confounders then add confounders
# Varying-effects for date, household, beach (hierarchical) 
# Use weakly informative priors for water contact

get_prior(agi3 ~ water_contact2 + (1 | beach/recruit_date/house_id), 
          family = bernoulli, data = data_follow)

priors <- c(set_prior("normal(0.3, 0.6)", class = "b", coef = "water_contact2Minimalcontact"),
            set_prior("normal(0.5, 0.5)", class = "b", coef = "water_contact2Bodyimmersion"),
            set_prior("normal(0.7, 0.4)", class = "b", coef = "water_contact2Swallowedwater"),
            set_prior("exponential(1)", class = "sd"))

m1 <- brm(agi3 ~ water_contact2 + (1 | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow, prior = priors, control = list(adapt_delta = 0.9),
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, 
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m1)
get_variables(m1)
plot(m1)
pp_check(m1, ndraws=100)
pp_check(m1, type = "stat", stat = "mean")
mcmc_acf(m1, pars = vars(contains("b_")), lags = 10)
mcmc_pairs(m1, pars = vars(contains("b_")), diag_fun = "den", off_diag_fun = "hex")

# Check if treating exposure as monotonic variable fits better than as indicator variable

data_follow <- data_follow |> 
  mutate(water_contact3 = factor(water_contact2, ordered = T, 
                                 levels = c("No contact", "Minimal contact", "Body immersion", "Swallowed water")))

get_prior(agi3 ~ mo(water_contact3) + (1 | beach/recruit_date/house_id), 
          family = bernoulli, data = data_follow)

# Dirichlet prior for the monotonic effect
# Prior expectation of dose-response relationship - first examination probability distribution of prior

dirichlet <- brms::rdirichlet(n = 1000, alpha = c(1, 2, 3)) |> 
  data.frame() |> 
  mutate(draw = 1:n()) 

dirichlet |> 
  pivot_longer(-draw, names_to = "level", values_to = "proportion") |> 
  group_by(level) |> 
  summarize(avg_prop = mean(proportion))

# Set priors and run model

priors2 <- c(set_prior("normal(0,1)",class= "b"),
             set_prior("exponential(1)", class = "sd"),
             set_prior("dirichlet(c(1, 2, 3))", class = "simo", coef = "mowater_contact31"))

m1.2 <- brm(agi3 ~ mo(water_contact3) + (1 | beach/recruit_date/house_id),
            family = bernoulli, data = data_follow, prior = priors2,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m1.2)
get_variables(m1.2)
plot(m1.2)
pp_check(m1.2, ndraws=100)
pp_check(m1.2, type = "stat", stat = "mean")
mcmc_acf(m1, pars = vars(contains("b_")), lags = 10)
mcmc_pairs(m1.2, pars = vars(contains("b_")), diag_fun = "den", off_diag_fun = "hex")
conditional_effects(m1.2)

loo(m1, m1.2)

# Monotonic version fits better

### Add E. coli variable with interaction (logged and mean centered/standardized version)

m2 <- brm(agi3 ~ mo(water_contact3)*log_e_coli_s + (1 | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m2)
get_variables(m2)
plot(m2)
pp_check(m2, ndraws=100)
pp_check(m2, type = "stat", stat = "mean")

conditional_effects(m2, effects = "log_e_coli_s:water_contact3")
conditional_effects(m2, effects = "water_contact3")

loo(m2)

# Compare to model with no interaction

m2.1 <- brm(agi3 ~ mo(water_contact3) + log_e_coli_s + (1 | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m2.1)
get_variables(m2.1)
plot(m2.1)
pp_check(m2.1, ndraws=100)
pp_check(m2.1, type = "stat", stat = "mean")

conditional_effects(m2.1, effects = "log_e_coli_s")
conditional_effects(m2.1, effects = "water_contact3")

loo(m2, m2.1)

# Interaction model has better fit, also fits with DAG and is part of research question/hypothesis

### Add minimal adjustment set of confounders
# Can add more iterations to improve estimation of house_id level, but trouble with this parameter

data_follow$ethnicity <- C(data_follow$ethnicity, contr.treatment, base=9)

m3 <- brm(agi3 ~ mo(water_contact3)*log_e_coli_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + (1 | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 5000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m3, robust = TRUE)
get_variables(m3)
plot(m3)
pp_check(m3, ndraws=100)
pp_check(m3, type = "stat", stat = "mean")
pairs(m3)

conditional_effects(m3, effects = "e_coli_s:water_contact3")
conditional_effects(m3, effects = "water_contact3") -> fit
fit$water_contact3

loo(m3)

# Given very low posterior predictions and high Pareto K values due to ~50% of households having n=1
# Re-run model without the household cluster and instead include indicator for household size >1

m4 <- brm(agi3 ~ mo(water_contact3)*log_e_coli_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (1 | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m4, robust = TRUE)
get_variables(m4)
plot(m4)
pp_check(m4, ndraws=100)
pp_check(m4, type = "stat", stat = "mean")

conditional_effects(m4, effects = "log_e_coli_s:water_contact3")
conditional_effects(m4, effects = "water_contact3") -> fit
fit$water_contact3

loo(m3, m4)

# Compare to model with site included as a third hierarchical varying effect

m4.1 <- brm(agi3 ~ mo(water_contact3)*log_e_coli_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (1 | site/beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m4.1, robust = TRUE)
get_variables(m4.1)
plot(m4.1)
pp_check(m4.1, ndraws=100)
pp_check(m4.1, type = "stat", stat = "mean")

conditional_effects(m4.1, effects = "log_e_coli_s:water_contact3")
conditional_effects(m4.1, effects = "water_contact3") -> fit
fit$water_contact3

loo(m4, m4.1)

# Model without site has better fit

# Compare to model with highest single sample E. coli 

m4.2 <- brm(agi3 ~ mo(water_contact3)*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group + 
              (1 | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors2,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m4.2, robust = TRUE)
get_variables(m4.2)
plot(m4.2)
pp_check(m4.2, ndraws=100)
pp_check(m4.2, type = "stat", stat = "mean")

conditional_effects(m4.2, effects = "log_e_coli_max_s:water_contact3")
conditional_effects(m4.2, effects = "water_contact3")

loo(m4, m4.2)

# Compare to model with qPCR enterococci instead of E. coli

m5 <- brm(agi3 ~ mo(water_contact3)*log_entero_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (1 | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m5, robust = TRUE)
get_variables(m5)
plot(m5)
pp_check(m5, ndraws=100)
pp_check(m5, type = "stat", stat = "mean")

conditional_effects(m5, effects = "log_entero_s:water_contact3")
conditional_effects(m5, effects = "water_contact3")

# Compare to qPCR enterococci highest single sample value

m5.1 <- brm(agi3 ~ mo(water_contact3)*log_entero_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (1 | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m5.1, robust = TRUE)
get_variables(m5.1)
plot(m5.1)
pp_check(m5.1, ndraws=100)
pp_check(m5.1, type = "stat", stat = "mean")

conditional_effects(m5.1, effects = "log_entero_max_s:water_contact3")
conditional_effects(m5.1, effects = "water_contact3")

loo(m5, m5.1)


# Compare to model with MST human sewage biomarker instead of E. coli

m6 <- brm(agi3 ~ mo(water_contact3)*log_mst_human_s + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group +
              (1 | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors2,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m6, robust = TRUE)
get_variables(m6)
plot(m6)
pp_check(m6, ndraws=100)
pp_check(m6, type = "stat", stat = "mean")

conditional_effects(m6, effects = "log_mst_human_s:water_contact3")
conditional_effects(m6, effects = "water_contact3")

## Reproduce other FIB models with same number of observations for LOO comparisons

data_follow_entero <- data_follow |> 
  drop_na(any_of(c("log_e_coli_s", "log_entero_s")))

m4.2_comp <- brm(agi3 ~ mo(water_contact3)*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
                 other_rec_act + beach_exp_food + sand_contact + household_group + 
                 (1 | beach/recruit_date),
               family = bernoulli, data = data_follow_entero, prior = priors2,
               iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
               backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

loo(m4.2_comp, m5.1)


### Negative control analysis - E. coli

data_neg_control <- data_follow |> filter(water_contact3 == "No contact")

priors_nc <- c(set_prior("normal(0, 1)", class = "b"), 
              set_prior("exponential(1)", class = "sd"))

m.nc <- brm(agi3 ~ log_e_coli_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (1 | beach/recruit_date),
          family = bernoulli, data = data_neg_control, prior = priors_nc,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m.nc, robust = TRUE)
plot(m.nc)
plot(m.nc, variable = "b_log_e_coli_s")
plot(m.nc, variable = "b_log_e_coli_s", combo = c("dens", "rank_overlay"))
pp_check(m.nc, ndraws=100)
pp_check(m.nc, type = "stat", stat = "mean")

conditional_effects(m.nc, effects = "log_e_coli_s")




