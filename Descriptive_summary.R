
pacman::p_load(
  rio,
  here,
  tidyverse,
  gtsummary,
  rstatix, 
  janitor, 
  flextable,
  viridis
)

# Load dataset

data <- import(here("data.xlsx"))
data_follow <- data |> filter(follow == "Yes") 

# Descriptive tables

data |> 
  select(age4, gender, education2, ethnicity, follow) |> 
  tbl_summary(by = follow, digits = list(all_categorical() ~ c(0, 1))) |> 
  add_overall() |>
  as_flex_table() 

data |> 
  select(cond_GI, other_rec_act, beach_exp_food, sand_contact, household_group, follow) |> 
  tbl_summary(by = follow, digits = list(all_categorical() ~ c(0, 1))) |> 
  add_overall() |>
  as_flex_table() 

data_follow |> tabyl(agi3)
data_follow |> tabyl(diarrhea3)

# FIB summary stats

old <- options(pillar.sigfig = 5)

data |> 
  distinct(recruit_date, .keep_all=TRUE) |> 
  get_summary_stats(e_coli, e_coli_max, entero_cce, entero_cce_max, turbidity) 

# Histograms

data |> group_by(date) |> ggplot(aes(x = e_coli)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = e_coli_s)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = log_e_coli)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = log_e_coli_s)) + geom_histogram()

data |> group_by(date) |> ggplot(aes(x = e_coli_max)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = entero_cce)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = entero_cce_max)) + geom_histogram()

data |> group_by(date) |> ggplot(aes(x = turbidity)) + geom_histogram()

# Examine FIB results by beach

data |>
  ggplot(aes(x = beach, y = e_coli, fill = beach)) +
  geom_violin() +
  geom_boxplot(width = 0.4, color="grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "E. coli geometric mean (CFU / 100 mL)", x = "Beach") + 
  theme(legend.position = "none")

data |>
  ggplot(aes(x = beach, y = e_coli_max, fill = beach)) +
  geom_violin() +
  geom_boxplot(width = 0.4, color="grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "E. coli highest single sample (CFU / 100 mL)", x = "Beach") + 
  theme(legend.position = "none")

data |>
  filter(site != "Toronto") |> 
  ggplot(aes(x = beach, y = entero_cce, fill = beach)) +
  geom_violin() +
  geom_boxplot(width = 0.4, color="grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "Enterococci geometric mean (CCE / 100 mL)", x = "Beach") + 
  theme(legend.position = "none")

data |>
  filter(site != "Toronto") |> 
  ggplot(aes(x = beach, y = entero_cce_max, fill = beach)) +
  geom_violin() +
  geom_boxplot(width = 0.4, color="grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "Enterococci highest single sample (CCE / 100 mL)", x = "Beach") + 
  theme(legend.position = "none")

# Examine water contact and AGI outcome relationship

data_follow |> 
  select(agi3, water_contact2) |> 
  tbl_summary(by = water_contact2, digits = list(all_categorical() ~ c(0, 1)),
              type = all_categorical() ~ "categorical")

data_follow |> 
  ggplot(aes(x = agi3, y = e_coli, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "E. coli geometric mean") +
  facet_grid(~water_contact2)

data_follow |> 
  ggplot(aes(x = agi3, y = e_coli_max, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "E. coli highest single sample") +
  facet_grid(~water_contact2)

data_follow |> 
  ggplot(aes(x = agi3, y = entero_cce, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "Enterococci geometric mean") +
  facet_grid(~water_contact2)

data_follow |> 
  ggplot(aes(x = agi3, y = entero_cce_max, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "Enterococci highest single sample") +
  facet_grid(~water_contact2)

data_follow |> 
  ggplot(aes(x = agi3, y = mst_human, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "HF183 human sewage marker (DNA / 100mL)") +
  facet_wrap(~water_contact2)

data_follow |> 
  filter(water_contact == "Yes") |> 
  ggplot(aes(x = agi3, y = mst_human, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "HF183 human sewage marker (DNA / 100mL)")

data_follow |> 
  ggplot(aes(x = agi3, y = mst_gull, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  facet_wrap(~water_contact2) +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "Seagull biomarker (DNA / 100mL)")

data_follow |> 
  ggplot(aes(x = agi3, y = turbidity, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "Water Turbidity" ) +
  facet_wrap(~water_contact2)

# Plot other water contact variables

data_follow |> 
  ggplot(aes(x = agi3, y = water_time, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "Time in the water (min)")

data_follow |> 
  ggplot(aes(x = water_time, y = e_coli, colour = agi3)) +
  geom_point() +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Time in the water (min)",
       y = "E. coli Geometric Mean")




