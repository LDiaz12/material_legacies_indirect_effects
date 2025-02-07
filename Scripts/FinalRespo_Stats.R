library(tidyverse)
library(here)

final_respo_norm <- read_csv(here("Data", "RespoFiles", "Final", "RespoR_Normalized_FinalRates.csv"))
metadata <- read_csv(here("Data", "MO24BEAST_Metadata.csv"))

## source to pH data using here function 
source(here("Scripts", "pH_slope_delta_TA.R"))

pH_summary <- pH_clean %>%
  group_by(TREATMENT, TANKID) %>%
  summarise(pH_mean = mean(pH_dailymean, na.rm = TRUE), 
            pH_range = mean(pH_range, na.rm = TRUE))

final_respo_meta_join <- final_respo_norm %>%
  right_join(metadata) %>%
  drop_na() %>%
  select(TANK_NUM, TREATMENT, GP, R, NP) %>%
  left_join(pH_summary %>%
              rename(TANK_NUM = TANKID))

final_respo_plot <- final_respo_meta_join %>%
  ggplot(aes(x=TREATMENT, y = GP)) + 
  geom_point()
final_respo_plot

final_respo_R <- final_respo_meta_join %>%
  ggplot(aes(x=TREATMENT, y = R)) +
  geom_point()
final_respo_R


final_respo_NP <- final_respo_meta_join %>%
  ggplot(aes(x=TREATMENT, y = NP)) +
  geom_point()
final_respo_NP

final_respo_GP_R <- final_respo_meta_join %>%
  ggplot(aes(x=TREATMENT, y = GP/R)) + 
  geom_point()
final_respo_GP_R

final_respo_meta_join %>%
  ggplot(aes(x=pH_mean, y = GP)) + 
  geom_point()

final_respo_meta_join %>%
  ggplot(aes(x=pH_range, y = GP)) +
  geom_point()
