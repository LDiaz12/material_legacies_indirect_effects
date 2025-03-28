library(tidyverse)
library(here)

## first read in metadata, and then all final calculation sheets ##
metadata <- read_csv(here("Data", "MO24BEAST_Metadata.csv"))
metadata <- metadata %>%
  filter(!TREATMENT == "Pre")

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

# TA means per treatment 
TA_data <- read_csv(here("Data", "Chemistry", "Cleaned_TA_Data_per_Treatment.csv"))

# final respo rates 
final_respo_rates <- read_csv(here("Data", "RespoFiles", "Final", "RespoR_Normalized_FinalRates.csv"))
final_respo_rates <- final_respo_rates %>%
  select(-c(DATE, RUN_NUM, umol.sec.corr, CHAMBER, Temp.C)) %>%
  drop_na()
final_respo_rates$CORAL_NUM <- as.numeric(final_respo_rates$CORAL_NUM)

## now use full_join to combine all bio parameters to the Metadata sheet ##
metadata_physio_full <- metadata %>% 
  full_join(endos) %>%
  full_join(chlorophyll) %>%
  full_join(afdw_metadata) %>%
  full_join(final_respo_rates) %>%
  filter(!TREATMENT == "Pre", 
         !chla.ug.cm2<0, 
         !chlc2.ug.cm2<0,
         !NP<0)
        
#write_csv(metadata_physio_full, here("Data", "MO24BEAST_physio_metadata.csv"))

chem_metadata <- read_csv(here("Data", "chem_metadata.csv"))

chem_full <- chem_metadata %>%
  full_join(DOC_tank_means) %>% 
  full_join(pH_data) %>%
  full_join(TA_data)

write_csv(chem_full, here("Data", "chem_full.csv"))

final_respo_rates$CORAL_NUM <- as.numeric(final_respo_rates$CORAL_NUM)

metadata_full_resp <- metadata_full %>%
  right_join(final_respo_rates)

#write_csv(metadata_full_resp, here("Data", "MO24BEAST_Metadata_FULL.csv"))



