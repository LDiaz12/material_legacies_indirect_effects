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

# endosymbiont final data sheet 
endos <- read_csv(here("Data", "Endosymbionts", "endo_data_calculated.csv"))
endos <- endos %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, endo_per_cm2)

# chlorophyll final data sheet 
chlorophyll <- read_csv(here("Data", "Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_chl_full_data.csv"))

# DOC, pH, TA, NEP, and NEC data per tank and treatment 
chem_summary_data <- read_csv(here("Data", "Chemistry", "chem_summary_data.csv"))

# final respo rates 
final_respo_data <- read_csv(here("Data", "RespoFiles", "full_respo_data.csv"))
final_respo_data <- final_respo_data %>%
  filter(!GENOTYPE == "BLANK")
final_respo_data$CORAL_NUM <- as.numeric(final_respo_data$CORAL_NUM)

final_respo_rates <- read_csv(here("Data", "RespoFiles", "Final", "RespoR_Normalized_FinalRates.csv"))
final_respo_rates <- final_respo_rates %>%
  select(-c(DATE, RUN_NUM, umol.sec.corr, CHAMBER, Temp.C)) %>%
  drop_na()

## now use full_join to combine all bio parameters to the Metadata sheet ##
metadata_physio_full <- metadata %>% 
  full_join(endos) %>%
  full_join(chlorophyll) %>%
  full_join(afdw_metadata) %>%
  full_join(final_respo_data)

metadata_physio_full <- metadata_physio_full %>%
  select(-c(Chlc_norm, Chl_total, chlc2.ug.cm2))
        
#write_csv(metadata_physio_full, here("Data", "MO24BEAST_physio_metadata.csv"))

physio_metadata <- read_csv(here("Data", "MO24BEAST_physio_metadata.csv"))

chem_metadata <- chem_summary_data

final_respo_rates$CORAL_NUM <- as.numeric(final_respo_rates$CORAL_NUM)

metadata_full <- physio_metadata %>%
  full_join(chem_metadata) %>%
  filter(!TREATMENT == "Pre")

#write_csv(metadata_full, here("Data", "MO24BEAST_Metadata_FULL.csv"))

# t-tests and effect size plot # 
metadata_full <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))
t.test(endo_per_cm2 ~ pH_mean, data = metadata_full)


# chla and endo regression # 
chl_endos <- lm(chla.ug.cm2 ~ endo_per_cm2, data = physio_metadata)
check_model(chl_endos)
chl_endo_plot <- physio_metadata %>%
  ggplot(aes(x=endo_per_cm2, y = chla.ug.cm2, color = TREATMENT)) + 
  geom_point()
chl_endo_plot

chl_plot <- physio_metadata %>%
  mutate(chl_endo = chla.ug.cm2/endo_per_cm2) %>%
  filter(chl_endo < 1) %>%
  ggplot(aes(x = TREATMENT, y = chl_endo)) + 
  geom_point(alpha = 0.5) + 
  stat_summary()
chl_plot

physio_metadata <- physio_metadata %>%
  mutate(chl_endo = chla.ug.cm2/endo_per_cm2)

chl_endo_model <- lmer(log(chl_endo) ~ TREATMENT + (1|GENOTYPE), data = physio_metadata %>%
                         filter(chl_endo < 1))
check_model(chl_endo_model)
anova(chl_endo_model)
summary(chl_endo_model)

chem_physio_data <- physio_metadata %>%
  full_join(chem_full)

chl_plot2 <- chem_physio %>%
  mutate(chl_endo = chla.ug.cm2/endo_per_cm2) %>%
  filter(mean_tissue_biomass < 0.00075) %>%
  ggplot(aes(x = mean_tissue_biomass, y = R)) +
  geom_point(aes(color = as.factor(TANK_NUM))) +
  #coord_trans(y = "log") + 
  geom_smooth(method = "lm")
chl_plot2


chl_endo_PC_model <- lmer(log(endo_per_cm2) ~ DOC_rangemean + (1|GENOTYPE), data = chem_physio %>%
                            mutate(chl_endo = chla.ug.cm2/endo_per_cm2))
anova(chl_endo_PC_model)


