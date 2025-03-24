
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

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = mean(e_coli_max_s, na.rm=TRUE), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

predictions(m4.2, re_formula = NA, by = "water_contact3", type = "response", newdata = nd)

pred <- predictions(m4.2, re_formula = NA, by = "water_contact3", type = "response", newdata = nd) |> posterior_draws()

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

avg_comparisons(m4.2, re_formula = NA, variables = "water_contact3", newdata = nd)

mfx <- comparisons(m4.2, re_formula = NA, variables = "water_contact3", by = "water_contact3", 
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "mean(Body immersion) - mean(No contact)" = "Body immersion",
                           "mean(Swallowed water) - mean(No contact)" = "Swallowed water",
                           "mean(Minimal contact) - mean(No contact)" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Water Contact Effect on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(-0.5, 100) -> Fig2A

# Risk ratio scale

avg_comparisons(m4.2, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnoravg", transform = "exp")

mfx <- comparisons(m4.2, re_formula = NA, comparison = "lnor", transform = "exp", 
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

avg_comparisons(m4.2, re_formula = NA, variables = "water_contact3", newdata = nd, by = "gender")

mfx <- comparisons(m4.2, re_formula = NA, variables = "water_contact3", by = "gender",
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
  labs(x = "Effect of Water Contact on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-0.5, 100) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)

# Age specific estimates 

avg_comparisons(m4.2, re_formula = NA, variables = "water_contact3", newdata = nd, by = "age4")

mfx <- comparisons(m4.2, re_formula = NA, variables = "water_contact3", by = "age4",
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
  labs(x = "Effect of Water Contact on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-0.5, 50) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)

# Predicted probabilities of E. coli, conditional on water contact level
# Sequence E. coli by range of logged, standardized and centered variable then back-transform

quantile(data_follow$log_e_coli_max, na.rm = TRUE)
quantile(data_follow$log_e_coli_max_s, na.rm = TRUE)

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = seq(-2.1, 1.9, by = 0.1), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m4.2, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(e_coli = exp(log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE))) 

ggplot(pred, aes(x = e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(100, 500, 1000, 1500, 2000))  -> Fig5a

ggplot(pred, aes(x = e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(100, 500, 1000, 1500, 2000)) +
  facet_wrap(~ water_contact3)   -> Fig6a

# Predictions on log scale

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

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)


avg_slopes(m4.2, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd)

avg_slopes(m4.2, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_contact3")


### Marginal effects of E. coli, conditional on water contact, at specific cut-points

quantile(data_follow$e_coli_max, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99), na.rm = TRUE)
quantile(data_follow$log_e_coli_max_s, probs = c(0.25, 0.5,  0.75, 0.9, 0.95, 0.99), na.rm = TRUE)

# Cut-points of 50th, 75th & 95th percentiles

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = quantile(log_e_coli_max_s, probs = c(0.50, 0.75, 0.95), na.rm = TRUE), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m4.2, re_formula = NA, by = c("log_e_coli_max_s", "water_contact3"), 
                    type = "response", newdata = nd) |> posterior_draws()

pred <- pred |> mutate(e_coli = exp(log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE))) |> 
  mutate(e_coli = round(e_coli, digits = 0)) |> 
  mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = water_contact3, fill = factor(log_e_coli_max_s))) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Predicted AGI Incident Risk per 1000 Beachgoers",
       y = "",
       fill = "") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ factor(e_coli)) +
  xlim(0, 200) 


### Beach-specific posterior probabilities and contrasts

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_e_coli_max_s = mean(e_coli_s, na.rm=TRUE), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No",
            beach = c("English Bay Beach", "Kitsilano Beach", "Grand Beach West", "Grand Beach East",
                      "Sunnyside", "Marie Curtis")) 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m4.2, re_formula = ~ (1 | beach), type = "response", newdata = nd) |> posterior_draws()

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

avg_comparisons(m4.2, re_formula = ~ (1 | beach),
                variables = "water_contact3", newdata = nd, by = "beach")

mfx <- comparisons(m4.2, re_formula = ~ (1 | beach), variables = "water_contact3", by = "beach",
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1))

ggplot(mfx, aes(x = draw, y = beach, fill = beach)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Effect of Water Contact on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-0.5, 50) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)


### Marginal effects for qPCR enterococci model ###

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_entero_max_s = mean(entero_max_s, na.rm=TRUE), 
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
  xlim(-0.5, 70) -> Fig2A

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
  xlim(0, 6) -> Fig2B

Fig2B <- Fig2B + scale_y_discrete(labels = NULL)
Fig2 <- Fig2A + Fig2B
Fig2 + plot_annotation(tag_levels = 'A')

# Predicted probabilities of qPCR Enterococcus relationship, conditional on water contact level

quantile(data_follow$entero_cce_max, na.rm = TRUE)
quantile(data_follow$log_entero_max_s, na.rm = TRUE)

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_entero_max_s = seq(-0.4, 6.6, by = 0.2), 
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
  theme(legend.position = "bottom")  -> Fig5b

Fig5b <- Fig5b + scale_y_discrete(labels = NULL)
Fig5 <- Fig5a + Fig5b
Fig5 + plot_annotation(tag_levels = 'A')


ggplot(pred, aes(x = entero, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Enterococci Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)   -> Fig6b

Fig6b <- Fig6b + scale_y_discrete(labels = NULL)
Fig6 <- Fig6a + Fig6b
Fig6 + plot_annotation(tag_levels = 'A')


avg_slopes(m5.1, re_formula = NA, variables = "log_entero_max_s", newdata = nd)

avg_slopes(m5.1, re_formula = NA, variables = "log_entero_max_s", newdata = nd, by = "water_contact3")




######################################################
### Marginal effects for MST human marker mt model ###

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            log_mst_human_mt_max_s = mean(entero_max_s, na.rm=TRUE), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

predictions(m7.1, re_formula = NA, by = "water_contact3", type = "response", newdata = nd)

pred <- predictions(m7.1, re_formula = NA, by = "water_contact3", type = "response", newdata = nd) |> posterior_draws()

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

avg_comparisons(m7.1, re_formula = NA, variables = "water_contact3", newdata = nd)

mfx <- comparisons(m7.1, re_formula = NA, variables = "water_contact3", by = "water_contact3", 
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
  xlim(-0.5, 70) -> Fig2A

# Odds ratio scale

avg_comparisons(m7.1, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnoravg", transform = "exp")

mfx <- comparisons(m7.1, re_formula = NA, type = "link", variables = "water_contact3",   
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
  xlim(0, 6) -> Fig2B

Fig2B <- Fig2B + scale_y_discrete(labels = NULL)
Fig2 <- Fig2A + Fig2B
Fig2 + plot_annotation(tag_levels = 'A')

# Predicted probabilities of MST relationship, conditional on water contact level

quantile(data_follow$mst_human_mt_max, na.rm = TRUE)
quantile(data_follow$log_mst_human_mt_max_s, na.rm = TRUE)

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            log_mst_human_mt_max_s = seq(-1.52, 1.68, by = 0.12), 
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
  theme(legend.position = "bottom")  -> Fig5b

Fig5b <- Fig5b + scale_y_discrete(labels = NULL)
Fig5 <- Fig5a + Fig5b
Fig5 + plot_annotation(tag_levels = 'A')


ggplot(pred, aes(x = mst_human_mt, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Human Mitochondrial DNA Highest Single Sample",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3)   -> Fig6b

Fig6b <- Fig6b + scale_y_discrete(labels = NULL)
Fig6 <- Fig6a + Fig6b
Fig6 + plot_annotation(tag_levels = 'A')


avg_slopes(m7.1, re_formula = NA, variables = "log_mst_human_mt_max_s", newdata = nd)

avg_slopes(m7.1, re_formula = NA, variables = "log_mst_human_mt_max_s", newdata = nd, by = "water_contact3")



### Marginal effects for sensitivity analysis models

# Time in water (min) model

quantile(data_follow$water_time, probs = c(0.25, 0.5, 0.9, 0.95), na.rm = TRUE)
quantile(data_follow$water_time_s, probs = c(0.25, 0.5, 0.9, 0.95), na.rm = TRUE)

nd <- data_follow |> 
  data_grid(water_time_s = seq(-0.9, 11.9, by = 0.5),
            log_e_coli_max_s = mean(e_coli_max_s, na.rm=TRUE), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

avg_slopes(m_watertime, re_formula = NA, variables = "water_time_s", newdata = nd)

nd <- data_follow |> 
  data_grid(water_time_s = quantile(water_time_s, probs = c(0.50, 0.95), na.rm = TRUE),
            log_e_coli_max_s = seq(-2.1, 1.9, by = 0.1), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No") 

avg_slopes(m_watertime, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_time_s")

nd <- data_follow |> 
  data_grid(water_time_s = seq(-0.9, 11.9, by = 0.5),
            log_e_coli_max_s = seq(-2.1, 1.9, by = 0.1), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = "No")

avg_slopes(m_watertime, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd)










