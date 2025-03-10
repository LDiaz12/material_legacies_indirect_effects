library(tidyverse)
library(here)

## code to combine all biological parameter data per individual coral ##
## first read in metadata, and then all final calculation sheets ##
metadata <- read_csv(here("Data", "MO24BEAST_Metadata.csv"))

# surface area final calculation
surface_area <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA_calculated.csv"))
surface_area <- surface_area %>%
  select(-c(CORALID, weight1_g, weight2_g, weight_of_wax_g, date))

# afdw final data sheet 
afdw_metadata <- read_csv(here("Data", "Data_Raw", "Growth", "coral_mean_biomass_calculated.csv"))

# endosymbiont final data sheet 
endos <- read_csv(here("Data", "Endosymbionts", "endo_data_calculated.csv"))
endos <- endos %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, endo_per_cm2)

# chlorophyll final data sheet 
chlorophyll <- read_csv(here("Data", "Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_chl_full_data.csv"))

# pH data per tank and treatment 

pH_data <- read_csv(here("Data", "Chemistry", "Cleaned_pH_Data_per_Treatment.csv"))

# DOC means per tank data 
DOC_tank_means <- read_csv(here("Data", "DOC", "DOC_tank_means_data.csv"))

# final respo rates 
final_respo_rates <- read_csv(here("Data", "RespoFiles", "Final", "RespoR_Normalized_FinalRates.csv"))

# TA means per treatment 
TA_data <- read_csv(here("Data", "Chemistry", "Cleaned_TA_Data_per_Treatment.csv"))


## now use right_join to combine all bio parameters to the Metadata sheet and continually add on to it ##
metadata_full <- metadata %>%
  right_join(surface_area)

metadata_full <- metadata_full %>%
  right_join(endos)

metadata_full <- metadata_full %>%
  right_join(chlorophyll, by = "CORAL_NUM", "TREATMENT")

metadata_full <- metadata_full %>% 
  select(-c(GENOTYPE.y, TREATMENT.y, TANK_NUM.y, CORAL_TANK_NUM.y, SA_cm_2.y)) %>%
  rename(GENOTYPE = GENOTYPE.x,
         TREATMENT = TREATMENT.x, 
         TANK_NUM = TANK_NUM.x, 
         CORAL_TANK_NUM = CORAL_TANK_NUM.x, 
         SA_cm_2 = SA_cm_2.x)

metadata_full$TANK_NUM <- as.factor(metadata_full$TANK_NUM)

metadata_full <- metadata_full %>% 
  right_join(DOC_tank_means) %>%
  drop_na()

#write_csv(metadata_full, here("Data", "MO24BEAST_Metadata_FULL.csv"))

metadata_full <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))

metadata_full_b <- metadata_full %>%
  right_join(final_respo_rates)

metadata_full_b <- metadata_full %>%
  right_join(pH_data)

metadata_full_b <- metadata_full_b %>%
  right_join(TA_data)

#write_csv(metadata_full_b, here("Data", "MO24BEAST_Metadata_FULL.csv"))
