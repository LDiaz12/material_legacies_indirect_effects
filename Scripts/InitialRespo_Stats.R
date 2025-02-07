library(tidyverse)
library(here)

initial_respo_norm <- read_csv(here("Data", "RespoFiles", "Initial", "RespoR_Normalized_InitialRates.csv"))
metadata <- read_csv(here("Data", "MO24BEAST_Metadata.csv"))

initial_respo_meta_join <- initial_respo_norm %>%
  right_join(metadata, by = c("CORAL_NUM", "GENOTYPE")) %>%
  select(-c(TANK_NUM, CORAL_TANK_NUM)) %>%
  drop_na()

