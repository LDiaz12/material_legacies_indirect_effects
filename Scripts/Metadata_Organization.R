library(tidyverse)
library(here)
library(car)

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
afdw <- afdw_metadata %>%
  drop_na()

# endosymbiont final data sheet 
endos <- read_csv(here("Data", "Endosymbionts", "endo_data_calculated.csv"))
endos <- endos %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, endo_per_cm2)

# chlorophyll final data sheet 
chlorophyll <- read_csv(here("Data", "Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_chl_full_data.csv"))
chl <- chlorophyll %>%
  drop_na()

# DOC, pH, TA, NEP, and NEC data per tank and treatment 
chem_summary_data <- read_csv(here("Data", "Chemistry", "chem_summary_data.csv"))
full_carb_chem_data <- read_csv(here("Data", "Chemistry", "Full_Carb_Chem_Data.csv"))

# final respo rates 
final_respo_data <- read_csv(here("Data", "RespoFiles", "full_respo_data.csv"))
final_respo_data <- final_respo_data %>%
  filter(!GENOTYPE == "BLANK")
final_respo_data$CORAL_NUM <- as.numeric(final_respo_data$CORAL_NUM)

final_respo_rates <- read_csv(here("Data", "RespoFiles", "Final", "RespoR_Normalized_FinalRates.csv"))
final_respo_rates <- final_respo_rates %>%
  select(-c(DATE, RUN_NUM, umol.sec.corr, CHAMBER)) %>%
  drop_na()

## now use full_join to combine all bio parameters to the Metadata sheet ##

chlorophyll$CORAL_NUM <- as.numeric(chlorophyll$CORAL_NUM)
final_respo_rates$CORAL_NUM <- as.numeric(final_respo_rates$CORAL_NUM)

metadata_physio_full <- metadata %>% 
  full_join(endos) %>%
  full_join(chl) %>%
  full_join(afdw) %>%
  full_join(final_respo_rates) %>%
  filter(TREATMENT != "Pre")

        
#write_csv(metadata_physio_full, here("Data", "MO24BEAST_physio_metadata.csv"))

physio_metadata <- read_csv(here("Data", "MO24BEAST_physio_metadata.csv"))

chem_metadata <- chem_summary_data

final_respo_rates$CORAL_NUM <- as.numeric(final_respo_rates$CORAL_NUM)

metadata_full <- metadata_physio_full %>%
  full_join(chem_metadata) 

write_csv(metadata_full, here("Data", "MO24BEAST_Metadata_FULL.csv"))

# write separate metadata files for DAY and NIGHT # 
# DAY #
chem_summary_data_DAY <- read_csv(here("Data", "Chemistry", "chem_summary_data_DAY.csv"))
metadata_full_DAY <- physio_metadata %>%
  full_join(chem_summary_data_DAY) %>%
  filter(!TREATMENT == "Pre")
#write_csv(metadata_full_DAY, here("Data", "metadata_full_DAY.csv"))

# NIGHT #
chem_summary_data_NIGHT <- read_csv(here("Data", "Chemistry", "chem_summary_data_NIGHT.csv"))
metadata_full_NIGHT <- physio_metadata %>%
  full_join(chem_summary_data_NIGHT) %>%
  filter(!TREATMENT == "Pre")
#write_csv(metadata_full_NIGHT, here("Data", "metadata_full_NIGHT.csv"))

