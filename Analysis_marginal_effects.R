
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
  rstan,
  viridis,
  bayestestR
)

### Conditional and marginal effects with 'marginaleffects' R package

## Examine predictions for E. coli model first
# Examine posterior predictions of water contact exposure (predicted probabilities)
# Predictions integrate over cluster-level variables (re_formula = NULL) 

avg_pred <- avg_predictions(m4.1, variables = "water_contact3", re_formula = NULL)
avg_pred

pred <- posterior_draws(avg_pred)
pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = water_contact3)) +
  stat_halfeye(fill = "#440154", slab_alpha = 0.5, .width = 0.95) +
  labs(x = "Predicted AGI Incident Risk per 1000 Beachgoers",
    y = "Level of Water Contact",
    subtitle = "Posterior Probability Distributions by Level of Water Contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(0,80)

# Examine marginal effects/contrast of water contact exposure effect - probability scale

avg_comp <- avg_comparisons(m4.1, re_formula = NULL, variables = "water_contact3")
avg_comp

mfx <- posterior_draws(avg_comp)
mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = contrast)) +
  stat_halfeye(fill = "#440154", slab_alpha = .5, .width = 0.95)  +
  annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, 
    fill = "gray70", alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "gray40", linewidth = 0.6) +
  labs(x = "AGI Incident Risk per 1000 Beachgoers", y="Level of Water Contact") +
  theme_minimal() +
  theme(legend.position = "none") + 
  xlim(-10, 60)


avg_comp <- avg_comp |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 


avg_comp <- avg_comp |> mutate(
   label_text = sprintf("%.0f [%.0f, %.0f]", 
    estimate * 1000, conf.low * 1000, conf.high * 1000)) 

ggplot(mfx, aes(x = draw, y = contrast)) +
  stat_halfeye(fill = "#440154", slab_alpha = .5, .width = 0.95)  +
  annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, 
           fill = "gray70", alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "gray40", linewidth = 0.6) +
  labs(x = "AGI Incident Risk per 1000 Beachgoers", y="Level of Water Contact") +
  theme_minimal() +
  theme(legend.position = "none") + 
  xlim(-10, 60) +
  geom_text(data = avg_comp,
    aes(x = 25, y = contrast, label = label_text),
    vjust = -2.5, size = 4,  fontface = "bold", color = "black") 

ggsave("Fig1.tif", width = 3, height = 3, scale = 1.8, units = "in", dpi = 300)

# Check proportion of posterior that is greater than 0 and other values

mfx |> group_by(contrast) |> 
  summarize(proportion_0 = mean(draw > 0),
            proportion_1 = mean(draw > 1),
            proportion_5 = mean(draw > 5),
            proportion_10 = mean(draw > 10),
            proportion_20 = mean(draw > 20))

# Calculate region of practical equivalence (ROPE)

mfx |> 
  group_by(contrast) |> 
  summarize(rope_result = list(rope(draw, range = c(-1, 1), ci = 0.95))) |> 
  tidyr::unnest(rope_result)

# Population-averaged (marginal) adjusted risk ratios

avg_comp <- avg_comparisons(m4.1, re_formula = NULL, variables = "water_contact3",
                comparison = "lnratioavg", transform = "exp")
avg_comp

mfx <- posterior_draws(avg_comp)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Body immersion) / mean(No contact))" = "Body immersion",
                           "ln(mean(Swallowed water) / mean(No contact))" = "Swallowed water",
                           "ln(mean(Minimal contact) / mean(No contact))" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = contrast)) +
  stat_halfeye(fill = "#440154", slab_alpha = .5, .width = 0.95)  +
  annotate("rect", xmin = 0.90, xmax = 1.10, ymin = -Inf, ymax = Inf, 
           fill = "gray70", alpha = 0.4) +
  geom_vline(xintercept = c(0.90, 1.10), linetype = "dashed", colour = "gray40", linewidth = 0.6) +
  labs(x = "Risk Ratio", y="Level of Water Contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(0,5) 


avg_comp <- avg_comp |> mutate(
  label_text = sprintf("%.2f [%.2f, %.2f]", 
                       estimate, conf.low, conf.high)) 

avg_comp <- avg_comp |> 
  mutate(contrast = recode(contrast, "ln(mean(Body immersion) / mean(No contact))" = "Body immersion",
                           "ln(mean(Swallowed water) / mean(No contact))" = "Swallowed water",
                           "ln(mean(Minimal contact) / mean(No contact))" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = contrast)) +
  stat_halfeye(fill = "#440154", slab_alpha = .5, .width = 0.95)  +
  annotate("rect", xmin = 0.90, xmax = 1.10, ymin = -Inf, ymax = Inf, 
           fill = "gray70", alpha = 0.4) +
  geom_vline(xintercept = c(0.90, 1.10), linetype = "dashed", colour = "gray40", linewidth = 0.6) +
  labs(x = "Risk Ratio", y="Level of Water Contact") +
  theme_minimal() +
  theme(legend.position = "none") + 
  xlim(0, 5) +
  geom_text(data = avg_comp,
            aes(x = 2.5, y = contrast, label = label_text),
            vjust = -2.5, size = 4,  fontface = "bold", color = "black")


ggsave("Fig2.tif", width = 3, height = 3, scale = 1.8, units = "in", dpi = 300)


mfx |> 
  group_by(contrast) |> 
  summarize(rope_result = list(rope(draw, range = c(0.90, 1.10), ci = 0.95))) |> 
  tidyr::unnest(rope_result)


# Gender specific estimates - RR scale
# First examine baseline risks by gender and Risk Difference

avg_predictions(m4.1, variables = "water_contact3", by = "gender", re_formula = NULL)
avg_comparisons(m4.1, variables = "water_contact3", by = "gender", re_formula = NULL)

avg_comp <- avg_comparisons(m4.1, re_formula = NULL, variables = "water_contact3", by = "gender",
                            comparison = "lnratioavg", transform = "exp")
avg_comp

mfx <- posterior_draws(avg_comp)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Body immersion) / mean(No contact))" = "Body immersion",
                           "ln(mean(Swallowed water) / mean(No contact))" = "Swallowed water",
                           "ln(mean(Minimal contact) / mean(No contact))" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1))  |> 
  mutate(gender = recode(gender, "man/boy" = "Man/boy", "woman/girl" = "Woman/girl",
                         "fluid/trans" = "Fluid/trans"))

ggplot(mfx, aes(x = draw, y = gender, fill = contrast)) +
  stat_halfeye( slab_alpha = .5)  +
  annotate("rect", xmin = 0.90, xmax = 1.10, ymin = -Inf, ymax = Inf, 
           fill = "gray70", alpha = 0.4) +
  geom_vline(xintercept = c(0.90, 1.10), linetype = "dashed", colour = "gray40", linewidth = 0.6) +
  labs(x = "Effect of Water Contact on AGI Incident Risk per 1000 Beachgoers", y = "Gender Identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(0, 5) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)


mfx |> 
  group_by(contrast, gender) |> 
  summarize(rope_result = list(rope(draw, range = c(0.90, 1.10), ci = 0.95))) |> 
  tidyr::unnest(rope_result)



# Age specific estimates - RR
# First examine baseline risks by gender and RD

avg_predictions(m4.1, variables = "water_contact3", by = "age4", re_formula = NULL)
avg_comparisons(m4.1, variables = "water_contact3", by = "age4", re_formula = NULL)

avg_comp <- avg_comparisons(m4.1, re_formula = NULL, variables = "water_contact3", by = "age4",
                            comparison = "lnratioavg", transform = "exp")
avg_comp

mfx <- posterior_draws(avg_comp)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Body immersion) / mean(No contact))" = "Body immersion",
                           "ln(mean(Swallowed water) / mean(No contact))" = "Swallowed water",
                           "ln(mean(Minimal contact) / mean(No contact))" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1))  |> 
  mutate(age4 = fct_relevel(age4, "0-4", "5-9", "10-14", "15-19", "20+"))

ggplot(mfx, aes(x = draw, y = age4, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  annotate("rect", xmin = 0.90, xmax = 1.10, ymin = -Inf, ymax = Inf, 
           fill = "gray70", alpha = 0.4) +
  geom_vline(xintercept = c(0.90, 1.10), linetype = "dashed", colour = "gray40", linewidth = 0.6) +
  labs(x = "Effect of Water Contact on AGI Incident Risk per 1000 Beachgoers", y = "Age Group") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(0, 5) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)

mfx |> 
  group_by(contrast, age4) |> 
  summarize(rope_result = list(rope(draw, range = c(0.90, 1.10), ci = 0.95))) |> 
  tidyr::unnest(rope_result)


# Predicted probabilities of E. coli, conditional on water contact level
# Sequence E. coli by range of logged, standardized and centered variable then back-transform

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  reframe(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

avg_pred <- avg_predictions(m4.1, type = "response", re_formula = NULL, variables = list(
                              log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2),
                              water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
avg_pred

pred <- posterior_draws(avg_pred)

pred <- pred  |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                "Body immersion", "Swallowed water")) 

pred <- pred |> 
  mutate(e_coli = exp(log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE))) 

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
  labs(x = "Log E. coli",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, by = 0.05)) +
    facet_wrap(~ water_contact3)  -> Fig_ecoli


# Predicted median and 95% CI values of E. coli cut-points stratified by water contact 

e_coli_predictions <- pred |> 
  group_by(water_contact3, e_coli) |> 
  summarize(median = median(draw),
            lower = quantile(draw, 0.025),
            upper = quantile(draw, 0.975))

# Average slope

avg_comparisons(m4.1, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"))

avg_comparisons(m4.1, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"), by = "water_contact3")


### Marginal effects of E. coli, conditional on water contact, at specific cut-points

# Cut-points of 25th, 50th, 75th & 95th percentiles

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  reframe(quantile = scales::percent(c(0.25, 0.5, 0.75, 0.95)),
            e_coli_max = quantile(e_coli_max, na.rm=TRUE, c(0.25, 0.5, 0.75, 0.95)),
            log_e_coli_max_s = quantile(log_e_coli_max_s, na.rm=TRUE, c(0.25, 0.5, 0.75, 0.95)))

list <- data |> distinct(recruit_date, .keep_all = TRUE) |>
  reframe(log_e_coli_max_s = quantile(log_e_coli_max_s, na.rm=TRUE, c(0.25, 0.5, 0.75, 0.95)))
list <- as.list(list)

avg_predictions(m4.1, re_formula = NULL, variables = "water_contact3",
                newdata = datagrid(log_e_coli_max_s = list$log_e_coli_max_s,
                                   grid_type = "counterfactual"), by = "log_e_coli_max_s")

avg_comp <- avg_comparisons(m4.1, re_formula = NULL, variables = "water_contact3",
                            newdata = datagrid(log_e_coli_max_s = list$log_e_coli_max_s,
                                               grid_type = "counterfactual"),
                                               by = "log_e_coli_max_s")
avg_comp

mfx <- posterior_draws(avg_comp)

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
  labs(x = "Water Contact Effect on AGI Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ factor(e_coli)) +
  xlim(-15, 100)

ggplot(mfx, aes(x = draw, y = factor(e_coli), fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Water Contact Effect on AGI Incident Risk per 1000 Beachgoers", y = "E. coli Percentile Value (CFU/100 mL)") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast) +
  xlim(-15, 100)


# Average comparisons - RR

avg_comparisons(m4.1, re_formula = NULL, variables = "water_contact3",
                newdata = datagrid(log_e_coli_max_s = list$log_e_coli_max_s,
                                   grid_type = "counterfactual"),
                by = "log_e_coli_max_s", comparison = "lnratioavg", transform = "exp")

# Check proportion of posterior that is greater than 0 and other values

mfx |> group_by(contrast, log_e_coli_max_s) |> 
  summarize(proportion_0 = mean(draw > 0),
            proportion_1 = mean(draw > 1),
            proportion_5 = mean(draw > 5),
            proportion_10 = mean(draw > 10),
            proportion_20 = mean(draw > 20))


### Site-specific posterior probabilities and contrasts - RR scale

avg_comp <- avg_comparisons(m4.1, re_formula = ~ (1 | site), variables = "water_contact3", by = "site",
                            comparison = "lnratioavg", transform = "exp")
avg_comp

mfx <- posterior_draws(avg_comp)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Body immersion) / mean(No contact))" = "Body immersion",
                           "ln(mean(Swallowed water) / mean(No contact))" = "Swallowed water",
                           "ln(mean(Minimal contact) / mean(No contact))" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = site, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratio", y = "Site",
       subtitle = "Posterior Probability Distributions", fill = "Level of Water contact vs. No Contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ contrast) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 8)



### Marginal effects for qPCR enterococci model ###

avg_pred <- avg_predictions(m5.1, variables = "water_contact3", re_formula = NULL)
avg_pred

pred <- posterior_draws(avg_pred)
pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = water_contact3)) +
  stat_halfeye(fill = "#440154", slab_alpha = 0.5) +
  labs(x = "Predicted AGI Incident Risk per 1000 Beachgoers",
       y = "Level of Water Contact",
       subtitle = "Posterior Probability Distributions by Level of Water Contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(0,80)

# Risk differences

avg_comp <- avg_comparisons(m5.1, re_formula = NULL, variables = "water_contact3")
avg_comp

mfx <- posterior_draws(avg_comp)
mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Body immersion - No contact" = "Body immersion",
                           "Swallowed water - No contact" = "Swallowed water",
                           "Minimal contact - No contact" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = contrast)) +
  stat_halfeye(fill = "#440154", slab_alpha = .5)  +
  annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, 
           fill = "gray70", alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "gray40", linewidth = 0.6) +
  labs(x = "AGI Incident Risk per 1000 Beachgoers", y="Level of Water Contact",
       subtitle = "Posterior Probability Distributions") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-25, 50)

mfx |> 
  group_by(contrast) |> 
  summarize(rope_result = list(rope(draw, range = c(-1, 1), ci = 0.95))) |> 
  tidyr::unnest(rope_result)


# Risk ratio scale

avg_comp <- avg_comparisons(m5.1, re_formula = NULL, variables = "water_contact3",
                            comparison = "lnratioavg", transform = "exp")
avg_comp

mfx <- posterior_draws(avg_comp)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Body immersion) / mean(No contact))" = "Body immersion",
                           "ln(mean(Swallowed water) / mean(No contact))" = "Swallowed water",
                           "ln(mean(Minimal contact) / mean(No contact))" = "Minimal contact")) |> 
  mutate(contrast = fct_relevel(contrast, "Body immersion", after = 1)) 

ggplot(mfx, aes(x = draw, y = contrast)) +
  stat_halfeye(fill = "#440154", slab_alpha = .5)  +
  annotate("rect", xmin = 0.90, xmax = 1.10, ymin = -Inf, ymax = Inf, 
           fill = "gray70", alpha = 0.4) +
  geom_vline(xintercept = c(0.90, 1.10), linetype = "dashed", colour = "gray40", linewidth = 0.6) +
  labs(x = "Risk Ratios", y="Level of Water Contact",
       subtitle = "Posterior Probability Distributions") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(0,5) 


# Predicted probabilities of qPCR Enterococcus relationship, conditional on water contact level

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  reframe(log_entero_max_s = range(log_entero_max_s, na.rm=TRUE))

avg_pred <- avg_predictions(m5.1, type = "response", re_formula = NULL, variables = list(
              log_entero_max_s = seq(-2.014249, 3.1953, by = 0.4), 
              water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
avg_pred

pred <- posterior_draws(avg_pred)

pred <- pred  |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                     "Body immersion", "Swallowed water")) 

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
  labs(x = "Log Enterococci",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, by = 0.05)) +
  facet_wrap(~ water_contact3)   -> Fig_entero


avg_comparisons(m5.1, re_formula = NULL, variables = list(log_entero_max_s = "iqr"))

avg_comparisons(m5.1, re_formula = NULL, variables = list(log_entero_max_s = "iqr"), by = "water_contact3")


### Marginal effects for MST human marker mt model ###

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  reframe(log_mst_human_mt_max_s = mean(log_mst_human_mt_max_s, na.rm=TRUE))
list <- as.list(list)

avg_predictions(m7.1, variables = "water_contact3", re_formula = NULL)

avg_comparisons(m7.1, re_formula = NULL, variables = "water_contact3")

avg_comparisons(m7.1, re_formula = NULL, variables = "water_contact3",
                comparison = "lnratioavg", transform = "exp")


data |> distinct(recruit_date, .keep_all = TRUE) |> 
  reframe(log_mst_human_mt_max_s = range(log_mst_human_mt_max_s, na.rm=TRUE))

avg_pred <- avg_predictions(m7.1, type = "response", re_formula = NULL, variables = list(
  log_mst_human_mt_max_s = seq(-1.695368, 1.567463, by = 0.4), 
  water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
avg_pred

pred <- posterior_draws(avg_pred)

pred <- pred  |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                     "Body immersion", "Swallowed water")) 

pred <- pred |> 
  mutate(mst_human_mt = exp(log_mst_human_mt_max_s*sd(data_follow$log_mst_human_mt_max, na.rm=TRUE) + mean(data_follow$log_mst_human_mt_max, na.rm=TRUE))) 

pred <- pred |> 
  mutate(log_mst_human_mt = log_mst_human_mt_max_s*sd(data_follow$log_mst_human_mt_max, na.rm=TRUE) + mean(data_follow$log_mst_human_mt_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_mst_human_mt, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Human Mitochondrial DNA marker",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, by = 0.05)) +
  facet_wrap(~ water_contact3)   -> Fig_human_mt


avg_comparisons(m7.1, re_formula = NULL, variables = list(log_mst_human_mt_max_s = "iqr"))

avg_comparisons(m7.1, re_formula = NULL, variables = list(log_mst_human_mt_max_s = "iqr"), by = "water_contact3")


### Marginal effects for MST human sewage biomarker model ###


list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  reframe(log_mst_human_max_s = mean(log_mst_human_max_s, na.rm=TRUE))
list <- as.list(list)

avg_predictions(m6.1, variables = "water_contact3", re_formula = NULL)

avg_comparisons(m6.1, re_formula = NULL, variables = "water_contact3")

avg_comparisons(m6.1, re_formula = NULL, variables = "water_contact3",
                comparison = "lnratioavg", transform = "exp")


data |> distinct(recruit_date, .keep_all = TRUE) |> 
  reframe(log_mst_human_max_s = range(log_mst_human_max_s, na.rm=TRUE))

avg_pred <- avg_predictions(m6.1, type = "response", re_formula = NULL, variables = list(
  log_mst_human_max_s = seq(-0.8906186, 2.3137543, by = 0.4), 
  water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
avg_pred

pred <- posterior_draws(avg_pred)

pred <- pred  |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                     "Body immersion", "Swallowed water")) 

pred <- pred |> 
  mutate(log_mst_human = log_mst_human_max_s*sd(data_follow$log_mst_human_max, na.rm=TRUE) + mean(data_follow$log_mst_human_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_mst_human, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Human HF183 DNA marker",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, by = 0.05)) +
  facet_wrap(~ water_contact3) -> Fig_human


avg_comparisons(m6.1, re_formula = NULL, variables = list(log_mst_human_max_s = "iqr"))

avg_comparisons(m6.1, re_formula = NULL, variables = list(log_mst_human_max_s = "iqr"), by = "water_contact3")


### MST seagull marker model

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  reframe(log_mst_gull_max_s = mean(log_mst_gull_max_s, na.rm=TRUE))
list <- as.list(list)

avg_predictions(m8.1, variables = "water_contact3", re_formula = NULL)

avg_comparisons(m8.1, re_formula = NULL, variables = "water_contact3")

avg_comparisons(m8.1, re_formula = NULL, variables = "water_contact3",
                comparison = "lnratioavg", transform = "exp")


data |> distinct(recruit_date, .keep_all = TRUE) |> 
  reframe(log_mst_gull_max_s = range(log_mst_gull_max_s, na.rm=TRUE))

avg_pred <- avg_predictions(m8.1, type = "response", re_formula = NULL, variables = list(
  log_mst_gull_max_s = seq(-2.718343, 1.832380, by = 0.4), 
  water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
avg_pred

pred <- posterior_draws(avg_pred)

pred <- pred  |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                     "Body immersion", "Swallowed water")) 

pred <- pred |> 
  mutate(log_mst_gull_max = log_mst_gull_max_s*sd(data_follow$log_mst_gull_max, na.rm=TRUE) + mean(data_follow$log_mst_gull_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_mst_gull_max, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Seagull Gull4 DNA marker",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, by = 0.05)) +
  facet_wrap(~ water_contact3)  -> Fig_gull


avg_comparisons(m8.1, re_formula = NULL, variables = list(log_mst_gull_max_s = "iqr"))

avg_comparisons(m8.1, re_formula = NULL, variables = list(log_mst_gull_max_s = "iqr"), by = "water_contact3")


### Turbidity model

avg_predictions(m9, variables = "water_contact3", re_formula = NULL)

avg_comparisons(m9, re_formula = NULL, variables = "water_contact3")

avg_comparisons(m9, re_formula = NULL, variables = "water_contact3",
                comparison = "lnratioavg", transform = "exp")

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  reframe(log_turbidity_s = range(log_turbidity_s, na.rm=TRUE))


avg_pred <- avg_predictions(m9, type = "response", re_formula = NULL, variables = list(
  log_turbidity_s = seq(-1.149963, 2.982043, by = 0.4), 
  water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
avg_pred

pred <- posterior_draws(avg_pred)

pred <- pred  |> mutate(water_contact3 = fct_relevel(water_contact3, "Minimal contact", 
                                                     "Body immersion", "Swallowed water")) 

pred <- pred |> 
  mutate(log_turbidity = log_turbidity_s*sd(data_follow$log_turbidity_s, na.rm=TRUE) + mean(data_follow$log_turbidity_s, na.rm=TRUE)) 

ggplot(pred, aes(x = log_turbidity, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log Turbidity",
       y = "Predicted Probability of AGI",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, by = 0.05)) +
  facet_wrap(~ water_contact3)  -> Fig_turbidity


avg_comparisons(m9, re_formula = NULL, variables = list(log_turbidity_s = "iqr"))

avg_comparisons(m9, re_formula = NULL, variables = list(log_turbidity_s = "iqr"), by = "water_contact3")


## Combine FIB plots together

Fig_human <- Fig_human + theme(legend.position = "none")
Fig_human_mt <- Fig_human_mt + theme(legend.position = "none")
Fig_gull <- Fig_gull + theme(legend.position = "none")
Fig_ecoli <- Fig_ecoli + theme(legend.position = "none")

Fig_FIB <- Fig_ecoli + Fig_human + Fig_human_mt + Fig_gull + Fig_entero + Fig_turbidity
Fig_FIB + plot_annotation(tag_levels = 'A') + plot_layout(ncol = 2)

ggsave("Fig3.tif", width = 4, height = 4, scale = 2, units = "in", dpi = 300)

remove(Fig_ecoli, Fig_human, Fig_human_mt, Fig_gull, Fig_entero, Fig_turbidity)


### Marginal effects for sensitivity analysis models

# Time in water (min) model

quantile(data_follow$water_time, na.rm = TRUE)
quantile(data_follow$water_time_s, na.rm = TRUE)

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

avg_comparisons(m_watertime, re_formula = NULL, variables = list(water_time_s = "iqr"))

avg_comparisons(m_watertime, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"))


avg_pred <- avg_predictions(m_watertime, type = "response", re_formula = NULL, variables = list(
  water_time_s = seq(-0.6142013, 9.7476305, by = 0.4)))

pred <- posterior_draws(avg_pred)

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

quantile(data_follow$water_time_s, probs = c(0.25, 0.50, 0.75, 0.95), na.rm = TRUE)

avg_pred <- avg_predictions(m_watertime, type = "response", re_formula = NULL, variables = list(
  water_time_s = c(0.6142013, -0.3033464, 0.2147452, 1.8726383),
  log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.4)))
avg_pred

pred <- posterior_draws(avg_pred)

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

avg_slopes(m_watertime, re_formula = NULL, variables = "log_e_coli_max_s", by = "water_time_s")


# Alternative outcomes, follow-up, and prior models

avg_comparisons(m_body, re_formula = NULL, variables = "water_exp_body")

avg_comparisons(m_body, re_formula = NULL,  variables = "water_exp_body",
                comparison = "lnratioavg", transform = "exp")

avg_comparisons(m_diar, re_formula = NULL, variables = "water_contact3")

avg_comparisons(m_diar, re_formula = NULL,  variables = "water_contact3",
                comparison = "lnratioavg", transform = "exp")

avg_comparisons(m_3day, re_formula = NULL,  variables = "water_contact3")

avg_comparisons(m_3day, re_formula = NULL, variables = "water_contact3",
                comparison = "lnratioavg", transform = "exp")

avg_comparisons(m_5day, re_formula = NULL, variables = "water_contact3")

avg_comparisons(m_5day, re_formula = NULL, variables = "water_contact3",
                comparison = "lnratioavg", transform = "exp")

avg_comparisons(m_weak, re_formula = NULL, variables = "water_contact3")

avg_comparisons(m_weak, re_formula = NULL, variables = "water_contact3",
                comparison = "lnratioavg", transform = "exp")

avg_comparisons(m4.1na, re_formula = NULL, variables = "water_contact3")

avg_comparisons(m4.1na, re_formula = NULL, variables = "water_contact3",
                comparison = "lnratioavg", transform = "exp")

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))


avg_comparisons(m_body, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"))

avg_comparisons(m_body, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"), by = "water_exp_body")

avg_comparisons(m_diar, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"))

avg_comparisons(m_diar, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"), by = "water_contact3")

avg_comparisons(m_3day, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"))

avg_comparisons(m_3day, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"), by = "water_contact3")

avg_comparisons(m_5day, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"))

avg_comparisons(m_5day, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"), by = "water_contact3")

avg_comparisons(m_weak, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"))

avg_comparisons(m_weak, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"), by = "water_contact3")

avg_comparisons(m4.1na, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"))

avg_comparisons(m4.1na, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"), by = "water_contact3")



avg_pred <- avg_predictions(m_body, type = "response", re_formula = NULL, variables = list(
  log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2),
  water_exp_body = c("No", "Yes")))
pred <- posterior_draws(avg_pred)

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
  facet_wrap(~ water_exp_body)


avg_pred <- avg_predictions(m_diar, type = "response", re_formula = NULL, variables = list(
  log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2),
  water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
pred <- posterior_draws(avg_pred)

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


avg_pred <- avg_predictions(m_3day, type = "response", re_formula = NULL, variables = list(
  log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2),
  water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
pred <- posterior_draws(avg_pred)

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


avg_pred <- avg_predictions(m_5day, type = "response", re_formula = NULL, variables = list(
  log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2),
  water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
pred <- posterior_draws(avg_pred)

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


avg_pred <- avg_predictions(m_weak, type = "response", re_formula = NULL, variables = list(
  log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2),
  water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
pred <- posterior_draws(avg_pred)

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


avg_pred <- avg_predictions(m4.1na, type = "response", re_formula = NULL, variables = list(
  log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2),
  water_contact3 = c("Minimal contact", "Body immersion", "Swallowed water")))
pred <- posterior_draws(avg_pred)

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
  reframe(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

avg_pred <- avg_predictions(m.nc, type = "response", re_formula = NULL, variables = list(
  log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2)))

pred <- posterior_draws(avg_pred)

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




