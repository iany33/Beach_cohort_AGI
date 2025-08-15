
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

get_prior(agi3 ~ 0 + Intercept + water_contact2 + (1 | beach/recruit_date/house_id), 
          family = bernoulli, data = data_follow)

priors <- c(set_prior("normal(0.3, 0.6)", class = "b", coef = "water_contact2Minimalcontact"),
            set_prior("normal(0.5, 0.5)", class = "b", coef = "water_contact2Bodyimmersion"),
            set_prior("normal(0.7, 0.4)", class = "b", coef = "water_contact2Swallowedwater"),
            set_prior("normal(0, 1)", class = "b", coef = "Intercept"),
            set_prior("exponential(1)", class = "sd"))

m1 <- brm(agi3 ~ 0 + Intercept + water_contact2 + (1 | beach/recruit_date/house_id),
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

get_prior(agi3 ~ 0 + Intercept + mo(water_contact3) + (1 | beach/recruit_date/house_id), 
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

m1.2 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3) + (1 | beach/recruit_date/house_id),
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

m2 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_e_coli_s + (1 | beach/recruit_date/house_id),
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

### Add minimal adjustment set of confounders

data_follow$ethnicity <- C(data_follow$ethnicity, contr.treatment, base=9)

m3 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_e_coli_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + (1 | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m3, robust = TRUE)
get_variables(m3)
plot(m3)
pp_check(m3, ndraws=100)
pp_check(m3, type = "stat", stat = "mean")
pairs(m3)

conditional_effects(m3, effects = "log_e_coli_s:water_contact3")
conditional_effects(m3, effects = "water_contact3") -> fit
fit$water_contact3

loo(m3)

# Examine PAV-adjusted calibration plot diagnostic

y <- data_follow |> mutate(agi4 = if_else(agi3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_e_coli_s", "age4", "gender", "education2", "ethnicity", "cond_GI", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$agi4
yrep <- posterior_predict(m3, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

# Re-run model without the household cluster and instead include indicator for household size >1

m4 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_e_coli_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (1 | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m4, robust = TRUE)
get_variables(m4)
plot(m4)
pp_check(m4, ndraws=100)
pp_check(m4, type = "stat", stat = "mean")

conditional_effects(m4, effects = "log_e_coli_s:water_contact3")
conditional_effects(m4, effects = "water_contact3") -> fit
fit$water_contact3

# Examine PAV-adjusted calibration plot diagnostic

y <- data_follow |> mutate(agi4 = if_else(agi3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_e_coli_s", "age4", "gender", "education2", "ethnicity", "cond_GI", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$agi4
yrep <- posterior_predict(m4, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])


loo(m3, m4)

# Using the PAV-adjusted calibration plots, the model with household clustering is mis-specified
# Model with household clustering is over-confident in predicting high prob cases and under-confident for low prob cases
# Moving forward will use household group instead of house id clustering

# Compare to model with highest single sample E. coli 

m4.1 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group +
              (1 | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors2,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m4.1, robust = TRUE)
get_variables(m4.1)
plot(m4.1)
pp_check(m4.1, ndraws=100)
pp_check(m4.1, type = "stat", stat = "mean")

conditional_effects(m4.1, effects = "log_e_coli_max_s:water_contact3")
conditional_effects(m4.1, effects = "water_contact3")

# Examine PAV-adjusted calibration plot diagnostic

y <- data_follow |> mutate(agi4 = if_else(agi3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_e_coli_max_s", "age4", "gender", "education2", "ethnicity", "cond_GI", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$agi4
yrep <- posterior_predict(m4.1, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

loo(m4, m4.1)

# Compare to model with random slopes for water contact, given cross-level interaction term

priors_rslopes <- c(set_prior("normal(0,1)", class = "b"),
                    set_prior("exponential(1)", class = "sd"),
                    set_prior("dirichlet(c(1, 2, 3))", class = "simo", coef = "mowater_contact31"),
                    set_prior("lkj(2)", class = "cor"))

m4.1r <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
               other_rec_act + beach_exp_food + sand_contact + household_group +
               (mo(water_contact3) | beach/recruit_date),
             family = bernoulli, data = data_follow, prior = priors_rslopes,
             iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
             backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m4.1r, robust = TRUE)
get_variables(m4.1r)
plot(m4.1r)
pp_check(m4.1r, ndraws=100)
pp_check(m4.1r, type = "stat", stat = "mean")
stancode(m4.1r)

conditional_effects(m4.1r, effects = "log_e_coli_max_s:water_contact3")
conditional_effects(m4.1r, effects = "water_contact3") -> fit
fit$water_contact3

# Examine PAV-adjusted calibration plot diagnostic

y <- data_follow |> mutate(agi4 = if_else(agi3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_e_coli_max_s", "age4", "gender", "education2", "ethnicity", "cond_GI", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$agi4
yrep <- posterior_predict(m4.1r, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

loo(m4.1, m4.1r)

# Despite loo comparison not showing much difference, random slope model shows slightly more conservative predictions
# Use random-slope approach to avoid over estimating/predicting main effects

# Compare again to mean E. coli model with random slopes

m4r <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_e_coli_s + age4 + gender + education2 + ethnicity + cond_GI + 
             other_rec_act + beach_exp_food + sand_contact + household_group + 
             (mo(water_contact3) | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors_rslopes,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m4r, robust = TRUE)
get_variables(m4r)
plot(m4r)
pp_check(m4r, ndraws=100)
pp_check(m4r, type = "stat", stat = "mean")

conditional_effects(m4r, effects = "log_e_coli_s:water_contact3")
conditional_effects(m4r, effects = "water_contact3") 

loo(m4r, m4.1r)

# Compare to model with qPCR enterococci instead of E. coli

m5 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_entero_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (mo(water_contact3) | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors_rslopes,
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

m5.1 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_entero_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group +
              (mo(water_contact3) | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors_rslopes,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m5.1, robust = TRUE)
get_variables(m5.1)
plot(m5.1)
pp_check(m5.1, ndraws=100)
pp_check(m5.1, type = "stat", stat = "mean")

conditional_effects(m5.1, effects = "log_entero_max_s:water_contact3")
conditional_effects(m5.1, effects = "water_contact3")

# Examine PAV-adjusted calibration plot diagnostic

y <- data_follow |> mutate(agi4 = if_else(agi3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_entero_max_s", "age4", "gender", "education2", "ethnicity", "cond_GI", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$agi4
yrep <- posterior_predict(m5.1, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])


loo(m5, m5.1)


# Compare to model with MST human sewage biomarker HF183 instead of E. coli

m6 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_mst_human_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (mo(water_contact3) | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors_rslopes,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m6, robust = TRUE)
get_variables(m6)
plot(m6)
pp_check(m6, ndraws=100)
pp_check(m6, type = "stat", stat = "mean")

conditional_effects(m6, effects = "log_mst_human_s:water_contact3")
conditional_effects(m6, effects = "water_contact3")


m6.1 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_mst_human_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group +
              (mo(water_contact3) | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors_rslopes,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m6.1, robust = TRUE)
get_variables(m6.1)
plot(m6.1)
pp_check(m6.1, ndraws=100)
pp_check(m6.1, type = "stat", stat = "mean")

conditional_effects(m6.1, effects = "log_mst_human_max_s:water_contact3")
conditional_effects(m6.1, effects = "water_contact3")

# Examine PAV-adjusted calibration plot diagnostic

y <- data_follow |> mutate(agi4 = if_else(agi3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_mst_human_max_s", "age4", "gender", "education2", "ethnicity", "cond_GI", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$agi4
yrep <- posterior_predict(m6.1, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

loo(m6, m6.1)

# Given many non-detect days with 0s, check ordinal version based on quartiles

m6.2 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*mo(mst_human_max_cut) + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group +
              (mo(water_contact3) | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors_rslopes,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m6.2, robust = TRUE)
get_variables(m6.2)
plot(m6.2)
pp_check(m6.2, ndraws=100)
pp_check(m6.2, type = "stat", stat = "mean")

conditional_effects(m6.2, effects = "mst_human_max_cut:water_contact3")
conditional_effects(m6.2, effects = "water_contact3")

# Same trend is noted, so use continuous version winch avoids loss of information

# Compare to model with MST Human mitochondrial DNA marker

m7 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_mst_human_mt_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (mo(water_contact3) | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors_rslopes,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m7, robust = TRUE)
get_variables(m7)
plot(m7)
pp_check(m7, ndraws=100)
pp_check(m7, type = "stat", stat = "mean")

conditional_effects(m7, effects = "log_mst_human_mt_s:water_contact3")
conditional_effects(m7, effects = "water_contact3")

m7.1 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_mst_human_mt_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group +
              (mo(water_contact3) | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors_rslopes,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m7.1, robust = TRUE)
get_variables(m7.1)
plot(m7.1)
pp_check(m7.1, ndraws=100)
pp_check(m7.1, type = "stat", stat = "mean")

conditional_effects(m7.1, effects = "log_mst_human_mt_max_s:water_contact3")
conditional_effects(m7.1, effects = "water_contact3")

loo(m7, m7.1)

# Examine PAV-adjusted calibration plot diagnostic

y <- data_follow |> mutate(agi4 = if_else(agi3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_mst_human_mt_max_s", "age4", "gender", "education2", "ethnicity", "cond_GI", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$agi4
yrep <- posterior_predict(m7.1, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])


m7.2 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*mo(mst_human_mt_max_cut) + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group +
              (mo(water_contact3) | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors_rslopes,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m7.2, robust = TRUE)
get_variables(m7.2)
plot(m7.2)
pp_check(m7.2, ndraws=100)
pp_check(m7.2, type = "stat", stat = "mean")

conditional_effects(m7.2, effects = "mst_human_mt_max_cut:water_contact3")
conditional_effects(m7.2, effects = "water_contact3")


### Compare to model with MST Seagull marker 

m8 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_mst_gull_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (mo(water_contact3) | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors_rslopes,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m8, robust = TRUE)
get_variables(m8)
plot(m8)
pp_check(m8, ndraws=100)
pp_check(m8, type = "stat", stat = "mean")

conditional_effects(m8, effects = "log_mst_gull_s:water_contact3")
conditional_effects(m8, effects = "water_contact3")

m8.1 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_mst_gull_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group +
              (mo(water_contact3) | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors_rslopes,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m8.1, robust = TRUE)
get_variables(m8.1)
plot(m8.1)
pp_check(m8.1, ndraws=100)
pp_check(m8.1, type = "stat", stat = "mean")

conditional_effects(m8.1, effects = "log_mst_gull_max_s:water_contact3")
conditional_effects(m8.1, effects = "water_contact3")

# Examine PAV-adjusted calibration plot diagnostic

y <- data_follow |> mutate(agi4 = if_else(agi3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_mst_gull_max_s", "age4", "gender", "education2", "ethnicity", "cond_GI", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$agi4
yrep <- posterior_predict(m8.1, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])


loo(m8, m8.1)

# Compare to turbidity model

m9 <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_turbidity_s + age4 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + household_group +
            (mo(water_contact3) | beach/recruit_date), 
            family = bernoulli, data = data_follow, prior = priors_rslopes,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m9, robust = TRUE)
get_variables(m9)
plot(m9)
pp_check(m9, ndraws=100)
pp_check(m9, type = "stat", stat = "mean")

conditional_effects(m9, effects = "log_turbidity_s:water_contact3")
conditional_effects(m9, effects = "water_contact3")


y <- data_follow |> mutate(agi4 = if_else(agi3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_turbidity_s", "age4", "gender", "education2", "ethnicity", "cond_GI", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$agi4
yrep <- posterior_predict(m9, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,]) 

## Reproduce other FIB models with same number of observations for LOO comparisons

data_follow_entero <- data_follow |> 
  drop_na(any_of(c("log_e_coli_s", "log_entero_s")))

m4.2_comp <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
                   other_rec_act + beach_exp_food + sand_contact + 
                   (mo(water_contact3) | beach/recruit_date) + (1 | beach:recruit_date:house_id),
               family = bernoulli, data = data_follow_entero, prior = priors_rslopes,
               iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
               backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

loo(m4.2_comp, m5.1)

loo(m4.2, m6.1, m7.1, m8.1)


data_follow_any <- data_follow |> 
  drop_na(any_of(c("log_e_coli_s", "log_entero_s", "log_mst_human_s", "log_mst_human_mt_s")))

m6.1_comp <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_mst_human_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
                   other_rec_act + beach_exp_food + sand_contact + 
                   (mo(water_contact3) | beach/recruit_date) + (1 | beach:recruit_date:house_id),
                family = bernoulli, data = data_follow_any, prior = priors_rslopes,
                iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
                backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

loo(m4.2_comp, m5.1, m6.1_comp)


### Negative control analysis - E. coli

data_neg_control <- data_follow |> filter(water_contact3 == "No contact")

priors_nc <- c(set_prior("normal(0, 1)", class = "b"), 
              set_prior("exponential(1)", class = "sd"))

m.nc <- brm(agi3 ~ 0 + Intercept + log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
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

conditional_effects(m.nc, effects = "log_e_coli_max_s")



### Sensitivity analyses ###

# Alternative exposure: time spent in water (min), mean centered and standardized

data_follow <- data_follow |> 
  mutate(water_time_s = (water_time - mean(water_time, na.rm = TRUE)) / sd(water_time, na.rm = TRUE))

data_follow |> group_by(date) |> ggplot(aes(x = water_time)) + geom_histogram()
data_follow |> group_by(date) |> ggplot(aes(x = water_time_s)) + geom_histogram()

priors3 <- c(set_prior("normal(0,1)",class= "b"),
             set_prior("exponential(1)", class = "sd"),
             set_prior("lkj(2)", class = "cor"))

m_watertime <- brm(agi3 ~ 0 + Intercept + water_time_s*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group +
              (mo(water_time_s) | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors3,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_watertime, robust = TRUE)
get_variables(m_watertime)
plot(m_watertime)
pp_check(m_watertime, ndraws=100)
pp_check(m_watertime, type = "stat", stat = "mean")

conditional_effects(m_watertime, effects = "log_e_coli_max_s:water_time_s")
conditional_effects(m_watertime, effects = "water_time_s")

# Alternative outcome: diarrhea

m_diar <- brm(diarrhea3 ~ 0 + Intercept + mo(water_contact3)*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
                other_rec_act + beach_exp_food + sand_contact + household_group +
                (mo(water_contact3) | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors_rslopes,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_diar, robust = TRUE)
get_variables(m_diar)
plot(m_diar)
pp_check(m_diar, ndraws=100)
pp_check(m_diar, type = "stat", stat = "mean")

conditional_effects(m_diar, effects = "log_e_coli_max_s:water_contact3")
conditional_effects(m_diar, effects = "water_contact3")

# Alternative follow-up: 3 days

m_3day <- brm(agi_3day ~ 0 + Intercept + mo(water_contact3)*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
                other_rec_act + beach_exp_food + sand_contact + household_group +
                (mo(water_contact3) | beach/recruit_date),
                family = bernoulli, data = data_follow, prior = priors_rslopes,
                iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
                backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_3day, robust = TRUE)
get_variables(m_3day)
plot(m_3day)
pp_check(m_3day, ndraws=100)
pp_check(m_3day, type = "stat", stat = "mean")

conditional_effects(m_3day, effects = "log_e_coli_max_s:water_contact3")
conditional_effects(m_3day, effects = "water_contact3")


# Alternative follow-up: 5 days

m_5day <- brm(agi_5day ~ 0 + Intercept + mo(water_contact3)*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
                other_rec_act + beach_exp_food + sand_contact + household_group +
                (mo(water_contact3) | beach/recruit_date),
              family = bernoulli, data = data_follow, prior = priors_rslopes,
              iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
              backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_5day, robust = TRUE)
get_variables(m_5day)
plot(m_5day)
pp_check(m_5day, ndraws=100)
pp_check(m_5day, type = "stat", stat = "mean")

conditional_effects(m_5day, effects = "log_e_coli_max_s:water_contact3")
conditional_effects(m_5day, effects = "water_contact3")


# Stronger priors: normal(0,0.5)

priors_strong <- c(set_prior("normal(0,0.5)", class = "b"),
                   set_prior("exponential(1)", class = "sd"),
                   set_prior("dirichlet(c(1, 3, 6))", class = "simo", coef = "mowater_contact31"),
                   set_prior("lkj(2)", class = "cor"))

m_strong <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
                  other_rec_act + beach_exp_food + sand_contact + household_group +
                  (mo(water_contact3) | beach/recruit_date),
              family = bernoulli, data = data_follow, prior = priors_strong,
              iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.99),
              backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_strong, robust = TRUE)
get_variables(m_strong)
plot(m_strong)
pp_check(m_strong, ndraws=100)
pp_check(m_strong, type = "stat", stat = "mean")

conditional_effects(m_strong, effects = "log_e_coli_max_s:water_contact3")
conditional_effects(m_strong, effects = "water_contact3")

# Weaker priors: normal(0,2)

priors_weak <- c(set_prior("normal(0,2)", class = "b"),
                 set_prior("exponential(1)", class = "sd"),
                 set_prior("dirichlet(c(1, 1, 1))", class = "simo", coef = "mowater_contact31"),
                 set_prior("lkj(2)", class = "cor"))

m_weak <- brm(agi3 ~ 0 + Intercept + mo(water_contact3)*log_e_coli_max_s + age4 + gender + education2 + ethnicity + cond_GI + 
                other_rec_act + beach_exp_food + sand_contact + household_group +
                (mo(water_contact3) | beach/recruit_date),
                family = bernoulli, data = data_follow, prior = priors_weak,
                iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
                backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_weak, robust = TRUE)
get_variables(m_weak)
plot(m_weak)
pp_check(m_weak, ndraws=100)
pp_check(m_weak, type = "stat", stat = "mean")

conditional_effects(m_weak, effects = "log_e_coli_max_s:water_contact3")
conditional_effects(m_weak, effects = "water_contact3")






