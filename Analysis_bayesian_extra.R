
### Bayesian analysis of beach cohort outcome data

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

# Check different priors for random effects variance parameters

priors <- tibble(half_t = extraDistr::rht(n = 10000, nu = 3, sigma = 2.5),
                 exp = rexp(10000, rate = 1),
                 gamma = rgamma(10000, rate = 2, shape = 2)) 

priors |> ggplot(aes(x = half_t)) + geom_density(colour = "blue") +
  geom_density(aes(x = exp)) + 
  geom_density(aes(x = gamma), colour = "red") + 
  xlim(0, 10) + theme_minimal()

priors1 <- c(set_prior("normal(0.3, 0.6)", class = "b", coef = "water_contact2Minimalcontact"),
            set_prior("normal(0.5, 0.5)", class = "b", coef = "water_contact2Bodyimmersion"),
            set_prior("normal(0.7, 0.4)", class = "b", coef = "water_contact2Swallowedwater"),
            set_prior("exponential(1)", class = "sd"))

m1.1 <- brm(agi3 ~ water_contact2 + (1 | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow, prior = priors1, control = list(adapt_delta = 0.9),
          iter = 3000, chains = 4, cores = 4, warmup = 1000, seed = 123, 
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

loo(m1, m1.1)

# Check different hierarchical structure - remove household effects

m2.1 <- brm(agi3 ~ mo(water_contact3)*e_coli_s + (1 | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

loo(m2, m2.1)

# Add confounder instead for whether individual was part of household with >1 participant

m3a <- brm(agi3 ~ mo(water_contact3)*e_coli_s + age1 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group + (1 | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors2,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))


# Check difference if also using monotonic effects for age and education

data_follow <- data_follow |> 
  mutate(age2 = factor(age1, ordered = T, 
                       levels = c("0-4", "5-9", "10-14", "15-19", "20-39", "40-59", "60+")),
         education3 = factor(education2, ordered = T, 
                             levels = c("high school or less", "college/trades", "bachelors", "post-graduate")))

get_prior(agi3 ~ mo(water_contact3) + mo(age2) + gender + mo(education3) + ethnicity + cond_GI + 
            prev_act1 + beach_exp_food + sand_contact + e_coli_s + (1 | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow)

priors3 <- c(set_prior("normal(0,1)",class= "b"),
             set_prior("dirichlet(c(1, 2, 3))", class = "simo", coef = "mowater_contact31"),
             set_prior("dirichlet(c(1,1,1,1,1,1))", class = "simo", coef = "moage21"),
             set_prior("dirichlet(c(1, 1, 1))", class = "simo", coef = "moeducation31"))

m2.1 <- brm(agi3 ~ mo(water_contact3) + mo(age2) + gender + mo(education3) + ethnicity + cond_GI + 
              prev_act1 + beach_exp_food + sand_contact + e_coli_s + (1 | beach/recruit_date/house_id),
            family = bernoulli, data = data_follow, prior = priors3,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m2.1)
get_variables(m2.1)
plot(m2.1)
pp_check(m2.1, ndraws=100)
pp_check(m2.1, type = "stat", stat = "mean")

conditional_effects(m2.1, effects = "e_coli_s")
conditional_effects(m2.1, effects = "water_contact3")

loo(m2, m2.1)

# Model with confounders as dummy variables has better fit and faster convergence

# Model different hierarchical effects

m3.1 <- brm(agi3 ~ mo(water_contact3)*e_coli_s + age1 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + (1 | beach/recruit_date),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

pp_check(m3.1, ndraws=100)

conditional_effects(m3.1, effects = "water_contact3")

m3.2 <- brm(agi3 ~ mo(water_contact3)*e_coli_s + age1 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + beach + (1 | recruit_date/house_id),
          family = bernoulli, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

pp_check(m3.2, ndraws=100)

conditional_effects(m3.2, effects = "water_contact3")

m3.3 <- brm(agi3 ~ mo(water_contact3)*e_coli_s + age1 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + (1 | beach) + 
              (1 | recruit_date) + (1 | house_id),
            family = bernoulli, data = data_follow, prior = priors2,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

pp_check(m3.3, ndraws=100)

conditional_effects(m3.3, effects = "water_contact3")

loo(m3, m3.1, m3.2, m3.3)



# Compare again to using simpler model with no monotonic effects for water contact

priors3 <- c(set_prior("normal(0.3, 0.6)", class = "b", coef = "water_contact2Minimalcontact"),
             set_prior("normal(0.5, 0.5)", class = "b", coef = "water_contact2Bodyimmersion"),
             set_prior("normal(0.7, 0.4)", class = "b", coef = "water_contact2Swallowedwater"),
             set_prior("normal(0, 1)", class = "b"),
             set_prior("exponential(1)", class = "sd"))

m4.1 <- brm(agi3 ~ water_contact2*e_coli_s + age4 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + household_group +
              (1 | beach/recruit_date),
            family = bernoulli, data = data_follow, prior = priors3,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

conditional_effects(m4.1, effects = "e_coli_s:water_contact2")
conditional_effects(m4.1, effects = "water_contact2")

loo(m4, m4.1)



# Model outcome as Poisson instead of Bernoulli given rare status
# Offset outcome as rate per 1000 beachgoers

data_follow <- data_follow |> mutate(
  agi4 = if_else(agi3 == "Yes", 1, 0))

get_prior(agi4 ~ mo(water_contact3) + age1 + gender + education2 + ethnicity + cond_GI + 
            prev_act1 + beach_exp_food + sand_contact + e_coli_s + 
            (1 | beach/recruit_date/house_id), 
          family = poisson, data = data_follow)

m2.2 <- brm(agi4 ~ mo(water_contact3) + age1 + gender + education2 + ethnicity + cond_GI + 
              prev_act1 + beach_exp_food + sand_contact + e_coli_s + 
              (1 | beach/recruit_date/house_id),
            family = poisson, data = data_follow, prior = priors2,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m2.2)
get_variables(m2.2)
plot(m2.2)
pp_check(m2.2, ndraws=100)
pp_check(m2.2, type = "stat", stat = "mean")

conditional_effects(m2.2, effects = "e_coli_s")
conditional_effects(m2.2, effects = "water_contact3")

loo(m2, m2.2)

# Offset by number of beachgoers recruited per day to get daily rate

data_follow <- data_follow |> group_by(date) |> 
  summarize(participants_day = n()) |> 
  left_join(data_follow)

get_prior(agi4 ~ mo(water_contact3) + age1 + gender + education2 + ethnicity + cond_GI + 
            prev_act1 + beach_exp_food + sand_contact + e_coli_s + 
            log(participants_day) + (1 | beach/recruit_date/house_id), 
          family = poisson, data = data_follow)

m2.3 <- brm(agi4 ~ mo(water_contact3) + age1 + gender + education2 + ethnicity + cond_GI + 
            prev_act1 + beach_exp_food + sand_contact + e_coli_s + 
             log(participants_day) + (1 | beach/recruit_date/house_id),
          family = poisson, data = data_follow, prior = priors2,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m2.3)
get_variables(m2.3)
plot(m2.3)
pp_check(m2.3, ndraws=100)
pp_check(m2.3, type = "stat", stat = "mean")

conditional_effects(m2.3, effects = "e_coli_s")
conditional_effects(m2.3, effects = "water_contact3")

loo(m2, m2.2, m2.3)

### Examine model with E. coli values modeled as non-linear spline (thin-plate)

get_prior(agi3 ~ mo(water_contact3) + s(e_coli_s) + (1 | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow)

m2.1 <- brm(agi3 ~ mo(water_contact3) + s(e_coli_s) + (1 | beach/recruit_date/house_id),
            family = bernoulli, data = data_follow, prior = priors2,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m2.1)
get_variables(m2.1)
plot(m2.1)
pp_check(m2.1, ndraws=100)
pp_check(m2.1, type = "stat", stat = "mean")
plot(conditional_smooths(m2.1), ask = FALSE)

conditional_effects(m2.1, effects = "e_coli_s")
conditional_effects(m2.1, effects = "water_contact3")

# Compare to model with binned E. coli indicator variable

m3.3 <- brm(agi3 ~ mo(water_contact3)*ecoli_thresholds + age1 + gender + education2 + ethnicity + cond_GI + 
              other_rec_act + beach_exp_food + sand_contact + (1 | beach/recruit_date/house_id),
            family = bernoulli, data = data_follow, prior = priors2,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m3.3)
get_variables(m3.3)
plot(m3.3)
pp_check(m3.3, ndraws=100)
pp_check(m3.3, type = "stat", stat = "mean")

conditional_effects(m3.3, effects = "ecoli_thresholds:water_contact3")
conditional_effects(m3.3, effects = "water_contact3")

loo(m2, m3, m3.2, m3.3)



### Examine "swimming" as exposure vs level of water contact

m1.3 <- brm(agi3 ~ water_act_swim + (1 | beach/recruit_date/house_id),
            family = bernoulli, data = data_follow, prior = priors,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.9),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m1.3)
get_variables(m1.3)
plot(m1.3)
pp_check(m1.3, ndraws=100)
pp_check(m1.3, type = "stat", stat = "mean")
mcmc_acf(m1.3, pars = vars(contains("b_")), lags = 10)
mcmc_pairs(m1.3, pars = vars(contains("b_")), diag_fun = "den", off_diag_fun = "hex")
conditional_effects(m1.3)

loo(m1, m1.2, m1.3)


# Compare to model with swallowed water only as the predictor to simplify

get_prior(agi3 ~ water_exp_mouth*e_coli_s + age1 + gender + education2 + ethnicity + cond_GI + 
            prev_act1 + beach_exp_food + sand_contact + (1 | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow)

priors3 <- c(set_prior("normal(0.7, 0.4)", class = "b", coef = "water_exp_mouthYes"),
             set_prior("normal(0, 1)",class= "b"))

m3.3 <- brm(agi3 ~ water_exp_mouth*e_coli_s + age1 + gender + education2 + ethnicity + cond_GI + 
              prev_act1 + beach_exp_food + sand_contact + (1 | beach/recruit_date/house_id),
            family = bernoulli, data = data_follow, prior = priors3,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m3.3)
get_variables(m3.3)
plot(m3.3)
pp_check(m3.3, ndraws=100)
pp_check(m3.3, type = "stat", stat = "mean")

conditional_effects(m3.3, effects = "e_coli_s:water_exp_mouth")

loo(m2, m3, m3.3)


# Model E. coli as spline

m3.1 <- brm(agi3 ~ s(e_coli_s, by = interaction(water_contact2)) + water_contact2 + age1 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + (1 | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow, prior = priors3,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m3.1)
get_variables(m3.1)
plot(m3.1)
pp_check(m3.1, ndraws=100)
pp_check(m3.1, type = "stat", stat = "mean")

conditional_effects(m3.1, effects = "e_coli_s:water_contact2")
conditional_effects(m3.1, effects = "water_contact2")

loo(m3, m3.1)

### Compare to model with varying slopes for E. coli

get_prior(agi3 ~ mo(water_contact3)*e_coli_s + age1 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + 
            (1 + e_coli_s | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow)

priors4 <- c(set_prior("normal(0,1)",class = "b"),
             set_prior("exponential(1)", class = "sd"),
             set_prior("dirichlet(c(1, 2, 3))", class = "simo", coef = "mowater_contact31"),
             set_prior("lkj(2)", class = "cor"))

m4 <- brm(agi3 ~ mo(water_contact3)*e_coli_s + age1 + gender + education2 + ethnicity + cond_GI + 
            other_rec_act + beach_exp_food + sand_contact + 
            (1 + e_coli_s | beach/recruit_date/house_id),
          family = bernoulli, data = data_follow, prior = priors4,
          iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 123, control = list(adapt_delta = 0.95),
          backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m4)
get_variables(m4)
plot(m4)
pp_check(m4, ndraws=100)
pp_check(m4, type = "stat", stat = "mean")

conditional_effects(m4, effects = "e_coli_s:water_contact3")
conditional_effects(m4, effects = "water_contact3")

loo(m3, m4)




### Examine marginal effects - water contact - using emmeans

get_variables(m3)

pred <- m3 |> epred_draws(re_formula = NA, newdata = tibble(
  water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
  age1 = "0-4", gender = "woman/girl", ethnicity = "White", e_coli_s = 0,
  education2 = "bachelors", cond_GI = "No", prev_act1 = "No", beach_exp_food = "No", sand_contact = "No"))

pred <- pred |> mutate(water_contact3 = factor(water_contact3)) |> 
                mutate(water_contact3 =   
                         fct_relevel(water_contact3, c("No contact", "Minimal contact", "Body immersion", "Swallowed water")))

pred <- pred |> mutate(draw = .epred*1000)

ggplot(pred, aes(x = draw, y = water_contact3, fill = water_contact3)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted AGI Incidence per 1000 Beachgoers", y = "Level of Water Contact",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  xlim(0, 1)

mfx <- pred |> ungroup() |> select(-.row) |>  
  pivot_wider(names_from = "water_contact3", values_from = ".epred") |> 
  clean_names() |> 
  mutate(swallowed_effect = swallowed_water - no_contact,
         immersion_effect = body_immersion - no_contact,
         minimal_effect = minimal_contact - no_contact) 

mfx |> median_hdi(swallowed_effect) 



pacman::p_load(emmeans)

mfx <- m3 |> emmeans(~ water_contact3, re_formula = NA, epred = TRUE, at = list(e_coli_s = 0, 
      age1 = "0-4", gender = "woman/girl", ethnicity = "White", education2 = "bachelors", 
      cond_GI = "No", prev_act1 = "No", beach_exp_food = "No", sand_contact = "No")) |> 
    contrast(method = "revpairwise") |> gather_emmeans_draws()
  
mfx |> median_hdi() 


### Marginal effects with marginaleffects R package
### Examine marginal effects - water contact exposure

avg_comparisons(m3, variables = "water_contact3", comparison = "lnoravg", transform = "exp")

data_follow |> tabyl(house_id) |> arrange(desc(percent)) |> head()
data_follow |> tabyl(recruit_date) |> arrange(desc(percent)) |> head()

nd <- data.frame(water_contact2 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
                 e_coli_s = rep(mean(data_follow$e_coli_s, na.rm = TRUE), times = 4),
                 age1 = rep("0-4", times = 4), 
                 gender = rep("woman/girl", times = 4), 
                 ethnicity = rep("White", times = 4), 
                 education2 = rep("bachelors", times = 4), 
                 cond_GI = rep("No", times = 4), 
                 prev_act1 = rep("No", times = 4), 
                 beach_exp_food = rep("No", times = 4), 
                 sand_contact = rep("No", times = 4),
                 beach = rep("Grand Beach West", times = 4), 
                 house_id = rep("31302953", times = 4), 
                 recruit_day = rep("MB_2024-07-23", times = 4))

predictions(m3, re_formula = NA, type = "response", 
            newdata = datagrid(model = m3, water_contact3 = unique))

pred <- predictions(m3, re_formula = NA, type = "response", 
                    newdata = datagrid(model = m3, water_contact3 = unique)) |>  posteriordraws()

pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = water_contact2, fill = water_contact2)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted AGI Incidence per 1000 Beachgoers", y = "Level of Water Contact",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  xlim(0, 30)

comparisons(m3, re_formula = NA, variables = "water_contact3",
            newdata = datagrid(model = m3)) 

mfx <- comparisons(m3, re_formula = NA, variables = "water_contact3", 
                   newdata = datagrid(model = m3)) |> posteriordraws()

mfx <- mfx |> mutate(draw = draw*1000)

ggplot(mfx, aes(x = draw, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Effect of Water Contact on AGI Incidence per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  xlim(-0.5, 1)

mfx <- comparisons(m2, re_formula = NA, type = "link", variables = "water_contact3", 
                   newdata = datagrid(model = m2)) |> posteriordraws()

ggplot(mfx, aes(x = exp(draw), fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Odds Ratios of Water Contact Level vs. No Contact", y = "") +
  theme_minimal() +
  theme(legend.position = "bottom")  


# Marginal effects - E. coli

quantile(data_follow$e_coli, na.rm = TRUE)
quantile(data_follow$e_coli_s, na.rm = TRUE)

pred <- predictions(m3, re_formula = NA, type = "response",
                    newdata = datagrid(model = m3, water_contact3 = unique,
                                       e_coli_s = seq(-0.6, 3.9, by = 0.3))) |> 
  posteriordraws()

pred <- pred |> mutate(e_coli = round(e_coli_s*sd(data_follow$e_coli, na.rm=TRUE) + mean(data_follow$e_coli, na.rm=TRUE)))

ggplot(pred, aes(x = e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "E. coli Geometric Mean",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(100, 200, 300, 500, 1000))

ggplot(pred, aes(x = e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "E. coli Geometric Mean",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(100, 200, 300, 500, 1000)) +
  facet_wrap(~ water_contact3) 

mfx <- slopes(m3, re_formula = NULL, type = "response", variable = "e_coli_s",
              newdata = datagrid(water_contact3 = unique,
                                 e_coli_s = c(0.04, 1.7))) |> 
  posteriordraws()

mfx <- mfx |> mutate(e_coli = round(e_coli_s*sd(data_follow$e_coli) + mean(data_follow$e_coli)))

ggplot(mfx, aes(x = draw, fill = factor(e_coli))) +
  stat_halfeye() +
  labs(x = "Marginal Effect of E. coli Geometric Mean Values",
       y = "Posterior Density",
       fill = "") +
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap(~ contact) 

slopes(m3, re_formula = NULL, type = "response", variable = "e_coli_s",
       newdata = datagrid(water_contact3 = unique,
                          e_coli_s = c(0.06, 1.59)))

# Effect of Beach

pred <- predictions(m3, type = "response", allow_new_levels = TRUE,
                    newdata = datagrid(model = m3, water_contact3 = unique,
                                       beach = unique)) |>  posteriordraws()

pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = beach, fill = water_contact3)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted AGI Incidence per 1000 Beachgoers", y = "Beach",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  xlim(0, 0.1) +
  facet_wrap(~water_contact3)




