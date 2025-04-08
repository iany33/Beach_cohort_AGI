
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
  viridis,
  gganimate
)

# Create animated plot - E. coli relationship with AGI stratified by water contact

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

pred <- predictions(m4.2r, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) 

pred <- pred |>  
  group_by(water_contact3, log_e_coli, drawid) |> 
  summarize(draw = mean(draw))

pred <- pred |> 
  mutate(drawid = as.numeric(drawid)) |> 
  filter(drawid <= 50) |> 
  mutate(drawid = as.factor(drawid))

p <- pred |> 
  ggplot(aes(x = log_e_coli, y = draw)) +
  geom_line(aes(y = draw, group = drawid), color = "#3b528b") +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of AGI") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact3) +
  transition_states(drawid, 0, 1) +
  shadow_mark(past = TRUE, future = TRUE, alpha = 1/20, color = "gray50")

p_gif <- 
  animate(p, fps = 2, width = 800, height = 400, res = 150, renderer = gifski_renderer()) 

anim_save("e_coli_predictions.gif", animation = p_gif)
