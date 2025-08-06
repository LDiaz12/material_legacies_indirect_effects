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
            mean_blastate = mean(`BLASTATE VOL (ML)`)) # mean AFDW and tissue biomass PER coral 

#afdw_initial <- afdw_sa3 %>%
 # filter(TREATMENT == "Pre") %>%
  #group_by(GENOTYPE, mean_tissue_biomass) %>%
  #select(-c(CORAL_NUM, TREATMENT, mean_AFDW)) %>%
  #rename(initial_biomass = mean_tissue_biomass)

#afdw_data_full <- afdw_sa3 %>%
#  full_join(afdw_initial) %>%
 # filter(TREATMENT != "Pre") %>%
  #filter(mean_tissue_biomass < 40) # two values wrong 

ggplot(afdw_data_full) +
  geom_point(aes(x = initial_biomass, y = mean_tissue_biomass, color = TREATMENT))

# AFDW represents the total tissue biomass of the coral

### Plots for AFDW ### 
# Creating a custom color palette for coloring by genotype with more distinct colors
# Don't use this for talks! Just for personal purposes 
c24 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4")

### make AFDW plotdata  ###
afdw_plotdata <- afdw_data_full %>%
  group_by(TREATMENT) %>%
  summarize(tissuebiomass_mean = mean(mean_tissue_biomass, na.rm = TRUE),
            tissuebiomass_se = sd(mean_tissue_biomass, na.rm = TRUE)/sqrt(n()))
afdw_plotdata


# first plot out mean tissue biomass per treatment
afdw_plot <- afdw_data_full %>%
  ggplot(aes(x = TREATMENT, y = mean_tissue_biomass, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Mean Coral Tissue Biomass" ~ (g ~ mL^-1 ~ cm^-2)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  geom_jitter(data = afdw_data_full, aes(x = TREATMENT, y = mean_tissue_biomass), alpha = 0.7) +   
  theme(axis.title = element_text(size = 12, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_line(color = "gray")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15)) + 
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1), 
               geom = "errorbar", color = "black", width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_text(data = afdw_plotdata, 
            aes(x = TREATMENT, y = tissuebiomass_mean, 
                label = paste0("", round(tissuebiomass_mean, 4))),
            vjust = -1, hjust = 1.8, color = "black", size = 4)
afdw_plot
ggsave(plot = afdw_plot, filename = here("Output", "afdw_plot.png"), width = 9, height = 6)

## AFDW stats ## 
afdw_model <- lmer(mean_tissue_biomass ~ TREATMENT + (1|GENOTYPE), data=afdw_data_full)
check_model(afdw_model) # assumptions test reveals one outlier (remove values > 0.00075)
summary(afdw_model)
anova(afdw_model) 

write_csv(afdw_data_full, here("Data", "Data_Raw", "Growth", "coral_mean_biomass_calculated.csv"))
