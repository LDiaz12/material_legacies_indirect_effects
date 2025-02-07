library(tidyverse)
library(here)


## code to combine all biological parameter data per individual coral ##
## first read in metadata, and then all final calculation sheets ##
metadata <- read_csv(here("Data", "MO24BEAST_Metadata.csv"))
metadata <- metadata %>%
  select(-c(SAMPLE_ID))

surface_area <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA_calculated.csv"))
surface_area <- surface_area %>%
  select(-c(CORALID, date))

# afdw <- read_csv(here("Data", "Data_Raw", "Growth", "MO24BEAST_AFDW.csv")) # need to add blastates and calculate average biomass per coral
# so that there's 1 measurement instead of 3

endos <- read_csv(here("Data", "Endosymbionts", "MOBEAST_FCM_dataframe.csv"))
endos <- endos %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, INCLUDE_VOL_UL, ENDO_COUNT)
endos$CORAL_NUM <- as.numeric(endos$CORAL_NUM) #NAs introduced by coercion bc of blanks. ignore warning


chlorophyll <- read_csv(here("Data", "Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_chl_data.csv"))

## now use right_join to combine all bio parameters to the Metadata sheet and continually add on to it ##
metadata <- metadata %>%
  right_join(surface_area)

metadata <- metadata %>%
  right_join(endos)

metadata <- metadata %>%
  right_join(chlorophyll)

write_csv(metadata, here("Data", "MO24BEAST_Metadata.csv"))



