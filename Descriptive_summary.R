
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

data <- data |> 
  mutate(water_contact2 = case_when(
    water_contact == "No" ~ "No contact",
    water_exp_mouth == "Yes" ~ "Swallowed water",
    water_exp_body == "Yes" ~ "Body immersion",
    TRUE ~ "Minimal contact")) |>
  mutate(water_contact2 = as.factor(water_contact2)) |> 
  mutate(water_contact2 = fct_relevel(water_contact2, "No contact", "Minimal contact")) 

data_follow <- data |> filter(follow == "Yes") 

# Create new ordinal exposure variable

data_follow <- data_follow |> 
  mutate(water_contact3 = factor(water_contact2, ordered = T, 
                                 levels = c("No contact", "Minimal contact", "Body immersion", "Swallowed water")))

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

data_follow |> group_by(house_id) |> 
  mutate(n = n(), house_size = case_when(
    n == 1 ~ 1, n == 2 ~ 2, n == 3 ~ 3,
    n == 4 ~ 4, n == 5 ~ 5, n == 6 ~ 6)) |> 
  tabyl(house_size)

# FIB summary stats

old <- options(pillar.sigfig = 5)

data |> 
  distinct(recruit_date, .keep_all=TRUE) |> 
  get_summary_stats(e_coli, e_coli_max, entero_cce, entero_cce_max, 
                    mst_human, mst_human_max, mst_human_mt, mst_human_mt_max) 

data |> 
  distinct(recruit_date, .keep_all=TRUE) |> 
  summarize(mst_human_max = quantile(mst_human_max, probs = c(0, .25, .50, .75, 1)))

data |> 
  distinct(recruit_date, .keep_all=TRUE) |> 
  summarize(mst_human_mt_max = quantile(mst_human_mt_max, probs = c(0, .25, .50, .75, 1)))

data <- data |> mutate(mst_human_max_cut = case_when(
  mst_human_max == 0 ~ "0",
  (mst_human_max > 0 & mst_human_max <=132) ~ "1-132",
  (mst_human_max > 132 & mst_human_max <=415) ~ "133-415",
  mst_human_max > 415 ~ ">415")) |> 
  mutate(mst_human_max_cut = factor(mst_human_max_cut, ordered = T, 
               levels = c("0", "1-132", "133-415", ">415")))
    
data <- data |> mutate(mst_human_mt_max_cut = case_when(
  mst_human_mt_max == 0 ~ "0",
  (mst_human_mt_max > 0 & mst_human_mt_max <=113) ~ "1-113",
  (mst_human_mt_max > 113 & mst_human_mt_max <=245) ~ "114-245",
  mst_human_mt_max > 245 ~ ">245"))  |> 
  mutate(mst_human_mt_max_cut = factor(mst_human_mt_max_cut, ordered = T,
               levels = c("0", "1-113", "114-245", ">245")))                       

data |> distinct(date, mst_human_yn, site) |> 
  tabyl(mst_human_yn, site) |> 
  adorn_percentages("col") |> 
  adorn_totals("row") |> 
  adorn_ns() 

data |> distinct(date, mst_human_mt_yn, site) |> 
  tabyl(mst_human_mt_yn, site) |> 
  adorn_percentages("col") |> 
  adorn_totals("row") |> 
  adorn_ns() 

data |> distinct(date, mst_goose_yn, site) |> 
  tabyl(mst_goose_yn, site) |> 
  adorn_percentages("col") |> 
  adorn_totals("row") |> 
  adorn_ns() 

data |> distinct(date, mst_human_yn, beach) |> 
  tabyl(mst_human_yn, beach) |> 
  adorn_percentages("col") |> 
  adorn_totals("row") |> 
  adorn_ns() 

data |> distinct(date, mst_human_mt_yn, beach) |> 
  tabyl(mst_human_mt_yn, beach) |> 
  adorn_percentages("col") |> 
  adorn_totals("row") |> 
  adorn_ns() 

data |> distinct(date, mst_goose_yn, beach) |> 
  tabyl(mst_goose_yn, beach) |> 
  adorn_percentages("col") |> 
  adorn_totals("row") |> 
  adorn_ns() 

# Histograms

data |> group_by(date) |> ggplot(aes(x = e_coli)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = e_coli_s)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = log_e_coli)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = log_e_coli_s)) + geom_histogram()

data |> group_by(date) |> ggplot(aes(x = e_coli_max)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = log_e_coli_max)) + geom_histogram()

data |> group_by(date) |> ggplot(aes(x = entero_cce)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = log_entero)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = log_entero_s)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = entero_cce_max)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = log_entero_max)) + geom_histogram()

data |> group_by(date) |> ggplot(aes(x = mst_human_max)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = log_mst_human_max_s)) + geom_histogram()

data |> group_by(date) |> ggplot(aes(x = mst_human_mt_max)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = log_mst_human_mt_max_s)) + geom_histogram()

data |> group_by(date) |> ggplot(aes(x = log_mst_gull_s)) + geom_histogram()

data |> group_by(date) |> ggplot(aes(x = turbidity)) + geom_histogram()

# Examine FIB results by beach

data |>
  ggplot(aes(x = beach, y = log_e_coli, fill = beach)) +
  geom_violin() +
  geom_boxplot(width = 0.4, color="grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "Log E. coli mean (CFU / 100 mL)", x = "Beach") + 
  theme(legend.position = "none")

data |> 
  distinct(recruit_date, e_coli_max, beach) |> 
  ggplot(aes(x = beach, y = e_coli_max, fill = beach)) +
  geom_violin() +
  geom_point() +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "E. coli highest single sample (CFU / 100 mL)", x = "Beach") + 
  theme(legend.position = "none") +
  geom_hline(yintercept = 235, linetype = "dashed") +
  annotate("text", x = 0.55, y = 300, label = "BAV")

data |>
  filter(site != "Toronto") |> 
  ggplot(aes(x = beach, y = log_entero, fill = beach)) +
  geom_violin() +
  geom_boxplot(width = 0.4, color="grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "Log enterococci mean (CCE / 100 mL)", x = "Beach") + 
  theme(legend.position = "none")

data |>
  filter(site != "Toronto") |> 
  distinct(recruit_date, entero_cce_max, beach) |> 
  ggplot(aes(x = beach, y = entero_cce_max, fill = beach)) +
  geom_violin() +
  geom_point() +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "Enterococci highest single sample (CCE / 100 mL)", x = "Beach") + 
  theme(legend.position = "none") +
  geom_hline(yintercept = 235, linetype = "dashed") +
  annotate("text", x = 0.55, y = 300, label = "BAV")

data |>
  ggplot(aes(x = beach, y = log_mst_human, fill = beach)) +
  geom_violin() +
  geom_boxplot(width = 0.4, color="grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "Log MST Human HF183 Biomarker (DNA copies / 100 mL)", x = "Beach") + 
  theme(legend.position = "none")

data |>
  ggplot(aes(x = beach, y = log_mst_human_mt, fill = beach)) +
  geom_violin() +
  geom_boxplot(width = 0.4, color="grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "Log MST Human Mitochondrial Biomarker (DNA copies / 100 mL)", x = "Beach") + 
  theme(legend.position = "none")

data |>
  ggplot(aes(x = beach, y = log_mst_gull, fill = beach)) +
  geom_violin() +
  geom_boxplot(width = 0.4, color="grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(y = "Log MST Gull4 Biomarker (DNA copies / 100 mL)", x = "Beach") + 
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
  ggplot(aes(x = agi3, y = log_e_coli_max, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "Log E. coli highest single sample") +
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
  ggplot(aes(x = agi3, y = log_mst_human, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "Log HF183 human sewage marker (DNA / 100mL)") +
  facet_wrap(~water_contact2)

data_follow |> 
  filter(water_contact == "Yes") |> 
  ggplot(aes(x = agi3, y = log_mst_human, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "HF183 human sewage marker (DNA / 100mL)")

data_follow |> 
  ggplot(aes(x = agi3, y = log_mst_human_mt, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "Log human mitochondrial marker (DNA / 100mL)") +
  facet_wrap(~water_contact2)

data_follow |> 
  ggplot(aes(x = agi3, y = log_mst_gull, fill = agi3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  facet_wrap(~water_contact2) +
  labs(x = "Acute gastrointestinal illness (AGI)",
       y = "Log seagull biomarker (DNA / 100mL)")

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




