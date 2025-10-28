
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
  cmdstanr,
  modelr,
  patchwork,
  viridis
)

### Conditional and marginal effects with 'marginaleffects' R package

## Examine predictions for E. coli model first
# Examine posterior predictions of water contact exposure (predicted probabilities)
# Predictions ignore cluster-level variables (re_formula = NA) to get overall averages

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

exp(list$log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = list$log_e_coli_max_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

predictions(m4.1r, re_formula = NA, by = "water_contact3", type = "response", newdata = nd)

pred <- predictions(m4.1r, re_formula = NA, by = "water_contact3", type = "response", newdata = nd) |> get_draws()

pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = water_contact3, fill = water_contact3)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted AGI Incident Risk per 1000 Beachgoers", y = "Level of Water Contact",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 150)  -> Fig1

# Examine marginal effects/contrast of water contact exposure effect - probability scale

avg_comparisons(m4.1r, re_formula = NA, variables = "water_contact3", newdata = nd)

mfx <- comparisons(m4.1r, re_formula = NA, variables = "water_contact3", by = "water_contact3", 
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Water Contact Effect on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(-15, 70) -> Fig2A

# Check proportion of posterior that is greater than 0

mfx |> group_by(contrast) |> 
  summarize(proportion_0 = mean(draw > 0),
            proportion_10 = mean(draw > 10))

# Population-averaged (marginal) adjusted risk ratios

avg_comparisons(m4.1r, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnratioavg", transform = "exp")

mfx <- comparisons(m4.1r, re_formula = NA, comparison = "lnratio", transform = "exp", 
                   variables = "water_contact3",   
                   by = "water_contact3", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Body immersion) / mean(No contact))" = "Body immersion",
                           "ln(mean(Swallowed water) / mean(No contact))" = "Swallowed water",
                           "ln(mean(Minimal contact) / mean(No contact))" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratios (vs. No Water Contact)", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0,6) -> Fig2B

Fig2B <- Fig2B + scale_y_discrete(labels = NULL)
Fig2 <- Fig2A + Fig2B
Fig2 + plot_annotation(tag_levels = 'A')

remove(Fig2, Fig2A, Fig2B)

# Gender specific estimates 

avg_comparisons(m4.1r, re_formula = NA, variables = "water_contact3", newdata = nd, by = "gender")

mfx <- comparisons(m4.1r, re_formula = NA, variables = "water_contact3", by = "gender",
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) |> 
  mutate(gender = recode(gender, "man/boy" = "Man/boy", "woman/girl" = "Woman/girl",
                         "fluid/trans" = "Fluid/trans"))

ggplot(mfx, aes(x = draw, y = gender, fill = gender)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect of Water Contact on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-10, 150) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)

# Age specific estimates 

avg_comparisons(m4.1r, re_formula = NA, variables = "water_contact3", newdata = nd, by = "age4")

mfx <- comparisons(m4.1r, re_formula = NA, variables = "water_contact3", by = "age4",
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) |> 
  mutate(age4 = fct_relevel(age4, "0-4", "5-9", "10-14", "15-19", "20+"))

ggplot(mfx, aes(x = draw, y = age4, fill = age4)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect of Water Contact on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-10, 150) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)

# Predicted probabilities of E. coli, conditional on water contact level
# Sequence E. coli by range of logged, standardized and centered variable then back-transform

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m4.1r, re_formula = NA, by = c("water_contact3", "log_e_coli_max_s"), 
                    type = "response", newdata = nd) |> get_draws()

pred <- pred |> 
  mutate(e_coli = exp(log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE))) 

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) 

ggplot(pred, aes(x = e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(100, 500, 1000, 1500, 2000)) 

ggplot(pred, aes(x = e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(100, 500, 1000, 1500, 2000)) +
  facet_wrap(~ water_contact3)

# Predictions on log scale

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom")

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
    facet_wrap(~ water_contact3)   -> Fig_ecoli

# Predicted median and 95% CI values of E. coli cut-points stratified by water contact 

e_coli_predictions <- pred |> 
  group_by(water_contact3, e_coli) |> 
  summarize(median = median(draw),
            lower = quantile(draw, 0.025),
            upper = quantile(draw, 0.975))

# Average slope

avg_comparisons(m4.1r, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd)

avg_comparisons(m4.1r, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd, by = "water_contact3")


### Marginal effects of E. coli, conditional on water contact, at specific cut-points

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(quantile = scales::percent(c(0.25, 0.5, 0.75, 0.95)),
            e_coli_max = quantile(e_coli_max, na.rm=TRUE, c(0.25, 0.5, 0.75, 0.95)),
            log_e_coli_max_s = quantile(log_e_coli_max_s, na.rm=TRUE, c(0.25, 0.5, 0.75, 0.95)))

list <- data |> distinct(recruit_date, .keep_all = TRUE) |>
  summarize(log_e_coli_max_s = quantile(log_e_coli_max_s, na.rm=TRUE, c(0.25, 0.5, 0.75, 0.95)))
list <- as.list(list)

# Cut-points of 25th, 50th, 75th & 95th percentiles

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = list$log_e_coli_max_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

mfx <- comparisons(m4.1r, re_formula = NA, variables = "water_contact3", by = "log_e_coli_max_s", 
                   newdata = nd) |> get_draws()

mfx <- mfx |> mutate(e_coli = exp(log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE))) |> 
  mutate(e_coli = round(e_coli, digits = 0)) |> 
  mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Minimal contact - No contact" = "Minimal contact",
                           "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = contrast, fill = factor(log_e_coli_max_s))) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Water Contact Effect on AGI Incident Risk per 1000 Beachgoers at Specific E. coli Values", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ factor(e_coli)) +
  xlim(-15, 150)

ggplot(mfx, aes(x = draw, y = factor(e_coli), fill = factor(log_e_coli_max_s))) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Water Contact Effect on AGI Incident Risk per 1000 Beachgoers at Specific E. coli Values", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast) +
  xlim(-15, 150)

# Average comparisons 

avg_comparisons(m4.1r, re_formula = NA, variables = "water_contact3", newdata = nd, by = "log_e_coli_max_s")

avg_comparisons(m4.1r, re_formula = NA, variables = "water_contact3", by = "log_e_coli_max_s",
                newdata = nd, comparison = "lnratioavg", transform = "exp")



### Site-specific posterior probabilities and contrasts

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = list$log_e_coli_max_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes",
            site = data_follow$site) 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m4.1r, re_formula = ~ (1 | site), 
                    type = "response", newdata = nd) |> get_draws()

pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = site, fill = site)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted AGI Incident Risk per 1000 Beachgoers", y = "Site",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ water_contact3) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 200) 

avg_comparisons(m4.1r, re_formula = ~ (1 | site),
                variables = "water_contact3", newdata = nd, by = "site")

mfx <- comparisons(m4.1r, re_formula = ~ (1 | site), variables = "water_contact3", by = "site",
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1))

ggplot(mfx, aes(x = draw, y = site, fill = site)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect of Water Contact on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-10, 100) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)


### Marginal effects for qPCR enterococci model ###

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_entero_max_s = mean(log_entero_max_s, na.rm=TRUE))

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_entero_max_s = mean(log_entero_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_entero_max_s = list$log_entero_max_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

predictions(m5.1, re_formula = NA, by = "water_contact3", type = "response", newdata = nd)

pred <- predictions(m5.1, re_formula = NA, by = "water_contact3", type = "response", newdata = nd) |> posterior_draws()

pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = water_contact3, fill = water_contact3)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted AGI Incident Risk per 1000 Beachgoers", y = "Level of Water Contact",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 250)

# Examine marginal effects/contrast of water contact exposure effect - probability scale

avg_comparisons(m5.1, re_formula = NA, variables = "water_contact3", newdata = nd)

mfx <- comparisons(m5.1, re_formula = NA, variables = "water_contact3", by = "water_contact3", 
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Water Contact Effect on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(-10, 100) 

# Risk ratio scale

avg_comparisons(m5.1, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnratioavg", transform = "exp")

mfx <- comparisons(m5.1, re_formula = NA, type = "link", variables = "water_contact3",   
                   by = "water_contact3", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1))

ggplot(mfx, aes(x = exp(draw), y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratios (vs. No Water Contact)", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 6) 


# Predicted probabilities of qPCR Enterococcus relationship, conditional on water contact level

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_entero_max_s = range(log_entero_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_entero_max_s = seq(-1.909019, 2.676555, by = 0.4), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m5.1, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(entero = exp(log_entero_max_s*sd(data_follow$log_entero_max, na.rm=TRUE) + mean(data_follow$log_entero_max, na.rm=TRUE))) 

pred <- pred |> 
  mutate(log_entero_max = log_entero_max_s*sd(data_follow$log_entero_max, na.rm=TRUE) + mean(data_follow$log_entero_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_entero_max, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Enterococci Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom")  

ggplot(pred, aes(x = log_entero_max, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Enterococci Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)   -> Fig_entero


avg_comparisons(m5.1, re_formula = NA, variables = list(log_entero_max_s = "iqr"), newdata = nd)

avg_comparisons(m5.1, re_formula = NA, variables = list(log_entero_max_s = "iqr"), newdata = nd, by = "water_contact3")


### Marginal effects for MST human marker mt model ###

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_mst_human_mt_max_s = range(log_mst_human_mt_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_mst_human_mt_max_s = seq(-1.356289, 1.741581, by = 0.4), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m7.1, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(mst_human_mt = exp(log_mst_human_mt_max_s*sd(data_follow$log_mst_human_mt_max, na.rm=TRUE) + mean(data_follow$log_mst_human_mt_max, na.rm=TRUE))) 

pred <- pred |> 
  mutate(log_mst_human_mt = log_mst_human_mt_max_s*sd(data_follow$log_mst_human_mt_max, na.rm=TRUE) + mean(data_follow$log_mst_human_mt_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_mst_human_mt, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Human Mitochondrial DNA Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)   -> Fig_human


avg_comparisons(m7.1, re_formula = NA, variables = list(log_mst_human_mt_max_s = "iqr"), newdata = nd)

avg_comparisons(m7.1, re_formula = NA, variables = list(log_mst_human_mt_max_s = "iqr"), newdata = nd, by = "water_contact3")


### Marginal effects for MST human sewage biomarker model ###

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_mst_human_max_s = range(log_mst_human_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_mst_human_max_s = seq(-1.179804, 1.995136, by = 0.4), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m6.1, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(log_mst_human = log_mst_human_max_s*sd(data_follow$log_mst_human_max, na.rm=TRUE) + mean(data_follow$log_mst_human_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_mst_human, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Human Sewage Biomarker Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom")

ggplot(pred, aes(x = log_mst_human, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Human Sewage Biomarker Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3) -> Fig_human_mt


avg_comparisons(m6.1, re_formula = NA, variables = list(log_mst_human_max_s = "iqr"), newdata = nd)

avg_comparisons(m6.1, re_formula = NA, variables = list(log_mst_human_max_s = "iqr"), newdata = nd, by = "water_contact3")


### MST seagull marker model

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_mst_gull_max_s = range(log_mst_gull_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_mst_gull_max_s = seq(-1.976815, 2.273201, by = 0.4), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m8.1, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(log_mst_gull_max = log_mst_gull_max_s*sd(data_follow$log_mst_gull_max, na.rm=TRUE) + mean(data_follow$log_mst_gull_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_mst_gull_max, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log MST Seagull Biomarker Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom")  

ggplot(pred, aes(x = log_mst_gull_max, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log MST Seagull Biomarker Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)  -> Fig_gull


avg_comparisons(m8.1, re_formula = NA, variables = list(log_mst_gull_max_S = "iqr"), newdata = nd)

avg_comparisons(m8.1, re_formula = NA, variables = list(log_mst_gull_max_s = "iqr"), newdata = nd, by = "water_contact3")


### Turbidity model

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_turbidity_s = range(log_turbidity_s, na.rm=TRUE))

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_turbidity_s = mean(log_turbidity_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_turbidity_s = list$log_turbidity_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

predictions(m9, re_formula = NA, by = "water_contact3", type = "response", newdata = nd)

avg_comparisons(m9, re_formula = NA, variables = "water_contact3", newdata = nd)

avg_comparisons(m9, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnratioavg", transform = "exp")


nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_turbidity_s = seq(-1.149963, 2.982043, by = 0.4), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m9, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(log_turbidity = log_turbidity_s*sd(data_follow$log_turbidity_s, na.rm=TRUE) + mean(data_follow$log_turbidity_s, na.rm=TRUE)) 

ggplot(pred, aes(x = log_turbidity, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Turbidity",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom")  

ggplot(pred, aes(x = log_turbidity, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Turbidity",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)  -> Fig_turbidity


avg_comparisons(m9, re_formula = NA, variables = list(log_turbidity_s = "iqr"), newdata = nd)

avg_comparisons(m9, re_formula = NA, variables = list(log_turbidity_s = "iqr"), newdata = nd, by = "water_contact3")


## Combine FIB plots together

Fig_human <- Fig_human + theme(legend.position = "none")
Fig_human_mt <- Fig_human_mt + theme(legend.position = "none")
Fig_gull <- Fig_gull + theme(legend.position = "none")
Fig_ecoli <- Fig_ecoli + theme(legend.position = "none")
Fig_turbidity <- Fig_turbidity + theme(legend.position = "none")

Fig5 <- Fig_ecoli + Fig_human + Fig_human_mt + Fig_gull + Fig_entero + Fig_turbidity
Fig5 + plot_annotation(tag_levels = 'A') + plot_layout(ncol = 2)

remove(Fig_ecoli, Fig_human, Fig_human_mt, Fig_gull, Fig_entero, Fig_turbidity, Fig5)



### Marginal effects for sensitivity analysis models

# Time in water (min) model

quantile(data_follow$water_time, na.rm = TRUE)
quantile(data_follow$water_time_s, na.rm = TRUE)

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_time_s = seq(-0.8186170, 8.5499957, by = 0.4),
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

avg_comparisons(m_watertime, re_formula = NA, variables = list(water_time_s = "iqr"), newdata = nd)


pred <- predictions(m_watertime, re_formula = NA, type = "response", newdata = nd) |> get_draws()

pred <- pred |> 
  mutate(water_time = water_time_s*sd(data_follow$water_time, na.rm=TRUE) + mean(data_follow$water_time, na.rm=TRUE))

ggplot(pred, aes(x = water_time, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Amount of Time in the Water (Min)",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") 


data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_time_s = quantile(water_time_s, probs = c(0.25, 0.50, 0.75, 0.95), na.rm = TRUE),
            log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

avg_slopes(m_watertime, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_time_s")


pred <- predictions(m_watertime, re_formula = NA, type = "response", newdata = nd) |> get_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) |> 
  mutate(water_time = water_time_s*sd(data_follow$water_time, na.rm=TRUE) + mean(data_follow$water_time, na.rm=TRUE))

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_time)


# Alternative outcomes, follow-up, and prior models

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = list$log_e_coli_max_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

avg_comparisons(m_diar, re_formula = NA, variables = "water_contact3", newdata = nd)

avg_comparisons(m_diar, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnratioavg", transform = "exp")

avg_comparisons(m_3day, re_formula = NA, variables = "water_contact3", newdata = nd)

avg_comparisons(m_3day, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnratioavg", transform = "exp")

avg_comparisons(m_5day, re_formula = NA, variables = "water_contact3", newdata = nd)

avg_comparisons(m_5day, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnratioavg", transform = "exp")

avg_comparisons(m_weak, re_formula = NA, variables = "water_contact3", newdata = nd)

avg_comparisons(m_weak, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnratioavg", transform = "exp")


data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

avg_comparisons(m_diar, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd)

avg_comparisons(m_diar, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd, by = "water_contact3")

avg_comparisons(m_3day, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd)

avg_comparisons(m_3day, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd, by = "water_contact3")

avg_comparisons(m_5day, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd)

avg_comparisons(m_5day, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), by = "water_contact3")

avg_comparisons(m_weak, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd)

avg_comparisons(m_weak, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd, by = "water_contact3")


pred <- predictions(m_diar, re_formula = NA, type = "response", newdata = nd) |> posterior_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of Diarrhea",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)

pred <- predictions(m_3day, re_formula = NA, type = "response", newdata = nd) |> posterior_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)


pred <- predictions(m_5day, re_formula = NA, type = "response", newdata = nd) |> posterior_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)

pred <- predictions(m_weak, re_formula = NA, type = "response", newdata = nd) |> posterior_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)


# Negative control

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes", 
            sand_contact = "No", household_group = "Yes") 

pred <- predictions(m.nc, re_formula = NA, type = "response", newdata = nd) |> posterior_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom")

