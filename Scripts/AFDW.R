### Script to read in AFDW data ###
### By Laurel Diaz ###
### Created on 10/21/24 #####

#Load libraries
library(tidyverse)
library(lubridate)
library(here)
library(janitor)
library(lme4)
library(lmerTest)
library(ggridges)
library(moments)
library(emmeans)
library(agricolae)
library(car)
library(performance)

### read in plate data ###
afdw_data <- read_csv(here("Data", "Data_Raw", "Growth", "MO24BEAST_AFDW.csv"))
afdw_data <- afdw_data %>%
  select(SAMPLEID, CORAL_NUM, GENOTYPE, `BLASTATE VOL (ML)`:AFDW)

metadata <- read_csv(here("Data", "MO24BEAST_Metadata.csv"))

afdw_data_1 <- full_join(afdw_data, metadata)

### read in surface area data ###
sa <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA_calculated.csv"))
### combine afdw and sa data sheets ### 
afdw_sa <- full_join(afdw_data_1, sa)

## New data frame with normalized afdw and mean AFDW per coral ##
afdw_sa2 <- afdw_sa %>%
  mutate(tissue_biomass = (AFDW*1000*`BLASTATE VOL (ML)`)/SA_cm_2) # normalize AFDW to surface area of the coral (g/ml/cm2 ) #


afdw_data_full <- afdw_sa2 %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, TANK_NUM, `BLASTATE VOL (ML)`, AFDW, SA_cm_2, tissue_biomass) %>%
  group_by(CORAL_NUM, GENOTYPE, TREATMENT) %>%
  summarize(mean_AFDW = mean(AFDW), # calculate means because of triplicates
            mean_tissue_biomass = mean(tissue_biomass),
            mean_blastate = mean(`BLASTATE VOL (ML)`)) %>%
  filter(TREATMENT != "Pre") %>%
  filter(mean_tissue_biomass < 40)

#write_csv(afdw_data_full, here("Data", "Data_Raw", "Growth", "coral_mean_biomass_calculated.csv"))
