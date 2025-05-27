
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

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = list$log_e_coli_max_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

predictions(m4.2r, re_formula = NA, by = "water_contact3", type = "response", newdata = nd)

pred <- predictions(m4.2r, re_formula = NA, by = "water_contact3", type = "response", newdata = nd) |> get_draws()

pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = water_contact3, fill = water_contact3)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted AGI Incident Risk per 1000 Beachgoers", y = "Level of Water Contact",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 150)

# Examine marginal effects/contrast of water contact exposure effect - probability scale

avg_comparisons(m4.2r, re_formula = NA, variables = "water_contact3", newdata = nd)

mfx <- comparisons(m4.2r, re_formula = NA, variables = "water_contact3", by = "water_contact3", 
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
  xlim(-10, 100) -> Fig2A

# Odds ratio scale

avg_comparisons(m4.2r, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnoravg", transform = "exp")

mfx <- comparisons(m4.2r, re_formula = NA, comparison = "lnor", transform = "exp", 
                   variables = "water_contact3",   
                   by = "water_contact3", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(odds(Body immersion) / odds(No contact))" = "Body immersion",
                           "ln(odds(Swallowed water) / odds(No contact))" = "Swallowed water",
                           "ln(odds(Minimal contact) / odds(No contact))" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Odds Ratios (vs. No Water Contact)", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0,8) -> Fig2B

Fig2B <- Fig2B + scale_y_discrete(labels = NULL)
Fig2 <- Fig2A + Fig2B
Fig2 + plot_annotation(tag_levels = 'A')

# Gender specific estimates 

avg_comparisons(m4.2r, re_formula = NA, variables = "water_contact3", newdata = nd, by = "gender")

mfx <- comparisons(m4.2r, re_formula = NA, variables = "water_contact3", by = "gender",
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
  xlim(-10, 100) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)

# Age specific estimates 

avg_comparisons(m4.2r, re_formula = NA, variables = "water_contact3", newdata = nd, by = "age4")

mfx <- comparisons(m4.2r, re_formula = NA, variables = "water_contact3", by = "age4",
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
  xlim(-10, 70) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)

# Predicted probabilities of E. coli, conditional on water contact level
# Sequence E. coli by range of logged, standardized and centered variable then back-transform

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = seq(-2.35034, 1.823747, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m4.2r, re_formula = NA, by = c("water_contact3", "log_e_coli_max_s"), 
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

avg_slopes(m4.2r, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd)

avg_slopes(m4.2r, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_contact3")


### Marginal effects of E. coli, conditional on water contact, at specific cut-points

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(quantile = scales::percent(c(0.25, 0.5, 0.75, 0.95)),
            e_coli_max = quantile(e_coli_max, c(0.25, 0.5, 0.75, 0.95)),
            log_e_coli_max_s = quantile(log_e_coli_max_s, c(0.25, 0.5, 0.75, 0.95)))

list <- data |> distinct(recruit_date, .keep_all = TRUE) |>
  summarize(log_e_coli_max_s = quantile(log_e_coli_max_s, c(0.25, 0.5, 0.75, 0.95)))
list <- as.list(list)

# Cut-points of 25th, 50th, 75th & 95th percentiles

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = list$log_e_coli_max_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

mfx <- comparisons(m4.2r, re_formula = NA, variables = "water_contact3", by = "log_e_coli_max_s", 
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(e_coli = exp(log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE))) |> 
  mutate(e_coli = round(e_coli, digits = 0)) |> 
  mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
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

avg_comparisons(m4.2r, re_formula = NA, variables = "water_contact3", newdata = nd, by = "log_e_coli_max_s")

avg_comparisons(m4.2r, re_formula = NA, variables = "water_contact3", by = "log_e_coli_max_s",
                newdata = nd, comparison = "lnoravg", transform = "exp")



### Beach-specific posterior probabilities and contrasts

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = list$log_e_coli_max_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No",
            beach = c("English Bay Beach", "Kitsilano Beach", "Grand Beach West", "Grand Beach East",
                      "Sunnyside", "Marie Curtis")) 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m4.2r, re_formula = ~ (1 | beach), type = "response", newdata = nd) |> posterior_draws()

pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = beach, fill = beach)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted AGI Incident Risk per 1000 Beachgoers", y = "Beach",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ water_contact3) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 300) 

avg_comparisons(m4.2r, re_formula = ~ (1 | beach),
                variables = "water_contact3", newdata = nd, by = "beach")

mfx <- comparisons(m4.2r, re_formula = ~ (1 | beach), variables = "water_contact3", by = "beach",
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1))

ggplot(mfx, aes(x = draw, y = beach, fill = beach)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect of Water Contact on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-10, 50) +
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
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

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
  xlim(0, 150)

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

# Odds ratio scale

avg_comparisons(m5.1, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnoravg", transform = "exp")

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
  labs(x = "Odds Ratios (vs. No Water Contact)", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 6) 


# Predicted probabilities of qPCR Enterococcus relationship, conditional on water contact level

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_entero_max_s = range(log_entero_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_entero_max_s = seq(-1.819538, 2.69703, by = 0.4), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m5.1, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(entero = exp(log_entero_max_s*sd(data_follow$log_entero_max, na.rm=TRUE) + mean(data_follow$log_entero_max, na.rm=TRUE))) 

ggplot(pred, aes(x = entero, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Enterococci Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") 

ggplot(pred, aes(x = entero, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Enterococci Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3) 

avg_slopes(m5.1, re_formula = NA, variables = "log_entero_max_s", newdata = nd)

avg_slopes(m5.1, re_formula = NA, variables = "log_entero_max_s", newdata = nd, by = "water_contact3")

# Log scale predictions

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



### Marginal effects for MST human marker mt model ###

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_mst_human_mt_max_s = range(log_mst_human_mt_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_mst_human_mt_max_s = seq(-1.3562894, 1.7415812, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m7.1, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(mst_human_mt = exp(log_mst_human_mt_max_s*sd(data_follow$log_mst_human_mt_max, na.rm=TRUE) + mean(data_follow$log_mst_human_mt_max, na.rm=TRUE))) 

ggplot(pred, aes(x = mst_human_mt, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Human Mitochondrial DNA Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)  

# Log scale predictions

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


avg_slopes(m7.1, re_formula = NA, variables = "log_mst_human_mt_max_s", newdata = nd)

avg_slopes(m7.1, re_formula = NA, variables = "log_mst_human_mt_max_s", newdata = nd, by = "water_contact3")


### Marginal effects for MST human sewage biomarker model ###

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_mst_human_max_s = range(log_mst_human_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_mst_human_max_s = seq(-1.179804, 1.995136, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

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

avg_slopes(m6.1, re_formula = NA, variables = "log_mst_human_max_s", newdata = nd)

avg_slopes(m6.1, re_formula = NA, variables = "log_mst_human_max_s", newdata = nd, by = "water_contact3")


### MST seagull marker model

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_mst_gull_max_s = range(log_mst_gull_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_mst_gull_max_s = seq(-1.976815, 2.273201, by = 0.4), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

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

avg_slopes(m8.1, re_formula = NA, variables = "log_entero_max_s", newdata = nd)

avg_slopes(m8.1, re_formula = NA, variables = "log_entero_max_s", newdata = nd, by = "water_contact3")


## Combine FIB plots together


Fig_human <- Fig_human + theme(legend.position = "none")
Fig_human_mt <- Fig_human_mt + theme(legend.position = "none")
Fig_gull <- Fig_gull + theme(legend.position = "none")
Fig_ecoli <- Fig_ecoli + theme(legend.position = "none")

Fig5 <- Fig_ecoli + Fig_human + Fig_human_mt + Fig_gull + Fig_entero
Fig5 + plot_annotation(tag_levels = 'A') + plot_layout(ncol = 2)


### Marginal effects for sensitivity analysis models

# Time in water (min) model

quantile(data_follow$water_time, na.rm = TRUE)
quantile(data_follow$water_time_s, na.rm = TRUE)

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(water_time_s = seq(-0.9, 11.9, by = 0.5),
            log_e_coli_max_s = list$log_e_coli_max_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

avg_slopes(m_watertime, re_formula = NA, variables = "water_time_s", newdata = nd)

nd <- data_follow |> 
  data_grid(water_time_s = quantile(water_time_s, probs = c(0.25, 0.50, 0.75, 0.95), na.rm = TRUE),
            log_e_coli_max_s = seq(-2.35034, 1.823747, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

avg_slopes(m_watertime, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_time_s")

nd <- data_follow |> 
  data_grid(water_time_s = seq(-0.9, 11.9, by = 0.5),
            log_e_coli_max_s = seq(-2.35034, 1.823747, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No")

avg_slopes(m_watertime, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd)

# Alternative outcomes, follow-up, and prior models

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = list$log_e_coli_max_s, 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

avg_comparisons(m_diar, re_formula = NA, variables = "water_contact3", newdata = nd)

avg_comparisons(m_diar, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnoravg", transform = "exp")

avg_comparisons(m_3day, re_formula = NA, variables = "water_contact3", newdata = nd)

avg_comparisons(m_3day, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnoravg", transform = "exp")

avg_comparisons(m_5day, re_formula = NA, variables = "water_contact3", newdata = nd)

avg_comparisons(m_5day, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnoravg", transform = "exp")

avg_comparisons(m_strong, re_formula = NA, variables = "water_contact3", newdata = nd)

avg_comparisons(m_strong, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnoravg", transform = "exp")

avg_comparisons(m_weak, re_formula = NA, variables = "water_contact3", newdata = nd)

avg_comparisons(m_weak, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnoravg", transform = "exp")


data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = seq(-2.35034, 1.823747, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

avg_slopes(m_diar, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd)

avg_slopes(m_diar, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_contact3")

avg_slopes(m_3day, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd)

avg_slopes(m_3day, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_contact3")

avg_slopes(m_5day, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd)

avg_slopes(m_5day, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_contact3")

avg_slopes(m_strong, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd)

avg_slopes(m_strong, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_contact3")

avg_slopes(m_weak, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd)

avg_slopes(m_weak, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_contact3")


pred <- predictions(m_diar, re_formula = NA, type = "response", newdata = nd) |> posterior_draws()

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


pred <- predictions(m_strong, re_formula = NA, type = "response", newdata = nd) |> posterior_draws()

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





