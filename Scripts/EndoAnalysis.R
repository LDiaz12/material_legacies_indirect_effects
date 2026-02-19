library(tidyverse)
library(here)
library(moments)
library(emmeans)
library(agricolae)

here()
endos <- read_csv(here("Data", "Endosymbionts", "MOBEAST_FCM_dataframe.csv"))
afdw_data <- read_csv(here("Data", "Data_Raw", "Growth", "MO24BEAST_AFDW.csv")) %>%
  select(CORAL_NUM, `BLASTATE VOL (ML)`) %>%
  group_by(CORAL_NUM) %>%
  summarize(mean_blastate = mean(`BLASTATE VOL (ML)`, na.rm = TRUE))


surface_area <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA_calculated.csv"))
metadata <- read_csv(here("Data", "MO24BEAST_Metadata.csv"))

## Endosymbiont counts are measured as event/uL

## Calculate mean endo counts per treatment ## 
endo_data <- endos %>%
  filter(TREATMENT != "BLANK") %>%
  select(CORAL_NUM, INCLUDE_VOL_UL, ENDO_COUNT)

## Calculate endo count per cm2 using blastate volume from "AFDW" data set ##
# combine endo_data with afdw data # 
endo_data$CORAL_NUM <- as.numeric(endo_data$CORAL_NUM)

endo_data_full <- endo_data %>%
  full_join(afdw_data %>%
              select(CORAL_NUM, mean_blastate)) %>%
  full_join(surface_area %>%
              select(CORAL_NUM, SA_cm_2)) %>%
  mutate(endo_per_cm2 = (ENDO_COUNT)*(1/1000)*(mean_blastate)*(1/SA_cm_2)) %>%
  full_join(metadata) %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, endo_per_cm2)

endo_data_full <- endo_data_full %>%
  filter(TREATMENT != "Pre")

#write_csv(endo_data_full, here("Data", "Endosymbionts", "endo_data_calculated.csv"))

