### Script to read in chlorophyll data from the Synergy 96 well plate in the Molecular Lab ####
### By Nyssa Silbiger ###
### Created on 6/14/2023 #####
### Edited by Laurel Diaz 9/12/24 ###

### load libraries ######
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
library(pals)
library(performance)
library(ggplot2)
#install.packages("ggpubr")
library(ggpubr)

### read in plate data ###
metadata <- read_csv(here("Data", "MO24BEAST_Metadata.csv"))

PlateData1<-read_csv(here("Data","Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_Chl_Run1_Plate1.csv"), skip = 39) #skips first 39 lines
PlateData2<-read_csv(here("Data","Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_Chl_Run1_Plate2.csv"), skip = 39)
PlateData2 <- PlateData2[-(24:96),] # deleting the rows with empty wells

MetaData1<-read_csv(here("Data","Data_Raw","Chl_Content", "Chl_Files", "Metadata1.csv")) 
MetaData2<-read_csv(here("Data","Data_Raw","Chl_Content", "Chl_Files", "Metadata2.csv"))
MetaData2 <- MetaData2 %>% 
  filter(!CORAL_NUM == "EMPTY") ## filter out empty wells 

Plate1 <- full_join(PlateData1, MetaData1)

Plate2 <- full_join(PlateData2, MetaData2)

## Combine both plates together ## 
plates_comb <- bind_rows(Plate1, Plate2)

sa <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA_calculated.csv"))


plate_comb_full <- plates_comb %>%
  select(-c(CORAL_TANK_NUM, TANK_NUM, TREATMENT, GENOTYPE)) %>%
  mutate(CORAL_NUM = as.numeric(CORAL_NUM)) %>%
  full_join(sa %>%
              select(CORAL_NUM, SA_cm_2)) %>%
  drop_na(CORAL_NUM, WELL)


##### Analysis #####

# Bring plate and metadata dataframes together #

full_data <- plate_comb_full %>% 
  mutate(adj_663 = `663`-`750`, # subtracting 750 is normalizing; absorbance value that 'should' be zero 
         adj_630 = `630` - `750`) %>%
  mutate(Chla_ug_ml = (11.43*adj_663)/0.6 - (0.64*adj_630)/0.6, # 0.6 is the path length 
         Chlc_ug_ml = (27.09*adj_630)/0.6 - (3.63*adj_663)/0.6) %>%
  mutate(chla_ug = Chla_ug_ml * BLASTATE_ML, # normalize to blastate volume
         chlc_ug = Chlc_ug_ml * BLASTATE_ML) %>%
  mutate(chla_ug_cm2 = chla_ug/SA_cm_2, # normalize to surface area 
         chlc_ug_cm2 = chlc_ug/SA_cm_2) %>%
  full_join(metadata)

### Chl a notes ###
# 11.43 is the extinction coefficient of chla at wavelength 663
# 0.64 is the extinction coefficient of chla at wavelength 630

### Chl c2 notes ###
# 27.09 is the extinction coefficient of chl c2 at wavelength 630
# 3.63 is the extinction coefficient of chl c2 at wavelength 663

### we subtract from wavelength of 750 because this is a correction for the turbidity of the sample
# 1ml sample
# 0.6 cm path length adjustment


## BELOW ONLY FOR METADATA FILE WRITING ##

chl_data_full <- full_data %>%
  filter(TREATMENT != "Pre") %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, chla_ug_cm2)

chl_full_filtered <- chl_data_full %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, chla_ug_cm2) %>%
  filter(chla_ug_cm2 < 7) %>% # three corals > 7 that are clear outliers 
  filter(chla_ug_cm2 > 0) # two negative corals 

#write_csv(chl_full_filtered, here("Data", "Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_chl_full_data.csv"))
