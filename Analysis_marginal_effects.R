
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
### Examine posterior predictions of water contact exposure (predicted probabilities)
# Predictions ignore cluster-level variables (re_formula = NA) to get overall averages

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            e_coli_s = mean(e_coli_s, na.rm=TRUE), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = 0) 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

predictions(m4, re_formula = NA, by = "water_contact3", type = "response", newdata = nd)

pred <- predictions(m4, re_formula = NA, by = "water_contact3", type = "response", newdata = nd) |> posterior_draws()

pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = water_contact3, fill = water_contact3)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted AGI Incident Risk per 1000 Beachgoers", y = "Level of Water Contact",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 200)

# Examine marginal effects/contrast of water contact exposure effect - probability scale

avg_comparisons(m4, re_formula = NA, variables = "water_contact3", newdata = nd)

mfx <- comparisons(m4, re_formula = NA, variables = "water_contact3", by = "water_contact3", 
                   newdata = nd) |> posteriordraws()

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

# Odds ratio scale

avg_comparisons(m4, re_formula = NA, variables = "water_contact3", newdata = nd,
                comparison = "lnoravg", transform = "exp")

mfx <- comparisons(m4, re_formula = NA, type = "link", variables = "water_contact3",   
                   by = "water_contact3", newdata = nd) |> posteriordraws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "mean(Body immersion) - mean(No contact)" = "Body immersion",
                           "mean(Swallowed water) - mean(No contact)" = "Swallowed water",
                           "mean(Minimal contact) - mean(No contact)" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = exp(draw), y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Odds Ratios of Water Contact Level vs. No Contact", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 10) -> Fig2B

Fig2B <- Fig2B + scale_y_discrete(labels = NULL)
Fig2 <- Fig2A + Fig2B
Fig2 + plot_annotation(tag_levels = 'A')

### Gender specific estimates 

avg_comparisons(m4, re_formula = NA, variables = "water_contact3", newdata = nd, by = "gender")

mfx <- comparisons(m4, re_formula = NA, variables = "water_contact3", by = "gender",
                   newdata = nd) |> posteriordraws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "mean(Body immersion) - mean(No contact)" = "Body immersion",
                           "mean(Swallowed water) - mean(No contact)" = "Swallowed water",
                           "mean(Minimal contact) - mean(No contact)" = "Minimal contact")) |> 
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

### Age specific estimates 

avg_comparisons(m4, re_formula = NA, variables = "water_contact3", newdata = nd, by = "age4")

mfx <- comparisons(m4, re_formula = NA, variables = "water_contact3", by = "age4",
                   newdata = nd) |> posteriordraws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "mean(Body immersion) - mean(No contact)" = "Body immersion",
                           "mean(Swallowed water) - mean(No contact)" = "Swallowed water",
                           "mean(Minimal contact) - mean(No contact)" = "Minimal contact")) |> 
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

### Predicted probabilities of E. coli, conditional on water contact level
# Sequence E. coli by range of standardized and centered variable then back-transform

quantile(data_follow$e_coli, na.rm = TRUE)
quantile(data_follow$e_coli_s, na.rm = TRUE)

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            e_coli_s = seq(-0.56, 4.40, by = 0.2), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = 0) 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m4, re_formula = NA, type = "response", newdata = nd) |> 
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
  scale_x_continuous(breaks = c(100, 500, 1000, 1500))

ggplot(pred, aes(x = e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "E. coli Geometric Mean",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(100, 500, 1000, 1500)) +
  facet_wrap(~ water_contact3) 


avg_slopes(m4, re_formula = NA, variables = "e_coli_s", newdata = nd)

avg_slopes(m4, re_formula = NA, variables = "e_coli_s", newdata = nd, by = "water_contact3")


### Marginal effects of E. coli, conditional on water contact, at specific cut-points

quantile(data_follow$e_coli, probs = c(0.5, 0.6, 0.75, 0.9, 0.95, 0.99), na.rm = TRUE)
quantile(data_follow$e_coli_s, probs = c(0.5, 0.6, 0.75, 0.9, 0.95, 0.99), na.rm = TRUE)

# Cut-points of 75th & 99th percentile

nd <- data_follow |> 
  data_grid(water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water"),
            e_coli_s = c(-0.08571, 3.213), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = 0) 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

mfx <- slopes(m4, re_formula = NA, type = "response", variable = "e_coli_s",
              newdata = nd) |> 
  posteriordraws()

mfx <- mfx |> mutate(e_coli = round(e_coli_s*sd(data_follow$e_coli, na.rm=TRUE) + mean(data_follow$e_coli, na.rm=TRUE)))

ggplot(mfx, aes(x = draw, y = water_contact3, fill = factor(e_coli))) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Marginal Effect of E. coli Geometric Mean Values",
       y = "Posterior Density",
       fill = "") +
  theme_classic() +
  theme(legend.position = "none") +
  xlim(-0.025, 0.05) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ factor(e_coli)) 


### Beach-specific posterior probabilities and contrasts

nd <- data_follow |> 
  data_grid(water_contact3 = c("No contact", "Minimal contact", "Body immersion", "Swallowed water"),
            e_coli_s = mean(e_coli_s, na.rm=TRUE), 
            age4 = c("0-4", "5-9", "10-14", "15-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity = "White", education2 = "bachelors", cond_GI = "No", other_rec_act = "Yes", 
            beach_exp_food = "Yes", sand_contact = "No", household_group = 0,
            beach = c("English Bay Beach", "Kitsilano Beach", "Grand Beach West", "Grand Beach East",
                      "Sunnyside", "Marie Curtis")) 

nd <- nd |> mutate(water_contact3 = fct_relevel(water_contact3, "No contact", "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- predictions(m4, re_formula = ~ (1 | beach), type = "response", newdata = nd) |> posterior_draws()

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

avg_comparisons(m4, re_formula = ~ (1 | beach),
                variables = "water_contact3", newdata = nd, by = "beach")

mfx <- comparisons(m4, re_formula = ~ (1 | beach), variables = "water_contact3", by = "beach",
                   newdata = nd) |> posteriordraws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "mean(Body immersion) - mean(No contact)" = "Body immersion",
                           "mean(Swallowed water) - mean(No contact)" = "Swallowed water",
                           "mean(Minimal contact) - mean(No contact)" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1))

ggplot(mfx, aes(x = draw, y = beach, fill = beach)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Effect of Water Contact on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-0.5, 50) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)



