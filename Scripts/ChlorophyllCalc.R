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
  mutate(CORAL_NUM = as.numeric(CORAL_NUM)) %>%
  full_join(sa) %>%
  drop_na(TREATMENT)


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
         chlc_ug_cm2 = chlc_ug/SA_cm_2)

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

chl_initial <- full_data %>%
  filter(TREATMENT == "Pre") %>%
  select(c(GENOTYPE, chla_ug_cm2)) %>%
  rename(initial_chla = chla_ug_cm2)

chl_data_full <- full_data %>%
  left_join(chl_initial)


chl_full_filtered <- chl_data_full %>%
  select(CORAL_NUM, TANK_NUM, TREATMENT, GENOTYPE, chla_ug_cm2, initial_chla) %>%
  filter(chla_ug_cm2 < 7) %>% # three corals > 7 that are clear outliers 
  filter(chla_ug_cm2 > 0) %>% # two negative corals 
  filter(TREATMENT != "Pre")

## PLOTS ## 
chl_full_filtered$TREATMENT <- factor(chl_full_filtered$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

chl_total_plot <- ggplot(chl_full_filtered, aes(x=TREATMENT, y=chla_ug_cm2, fill=TREATMENT)) + 
  geom_boxplot() + 
  labs(x="Treatment", y="Total Chl (µg/cm2)") + 
  scale_fill_manual(values=c("Control"="blue","Algae_Dom"="darkgreen",
                             "Coral_Dom" = "coral","Rubble_Dom" = "tan"), name="Treatment") +
  theme_bw(base_size=14) +
  geom_point(position = position_dodge(width=0.4))
chl_total_plot 


#ggsave(plot = chl_total_plot, filename = here("Output", "chl_total_plot.png"), width = 9, height = 6)
## PRE treatment signifies chl measurements from corals BEFORE experiment, so they were not put into a treatment

# Creating a custom color palette for coloring by genotype with more distinct colors
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


# chl a plot  #
chl_a_plot <- chl_full_filtered %>%
  ggplot(aes(x = TREATMENT, y = chla_ug_cm2, color = TREATMENT)) + # color by treatment for talks
  labs(x = "Treatment", y = expression(bold("Chlorophyll a" ~ (µg ~ cm^-2)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  geom_jitter(width = 0.1, alpha = 0.7) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  theme(legend.position = "none") +
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1), 
               geom = "errorbar", color = "black", width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                     "Rubble_Dom" = "tan"))
chl_a_plot
#ggsave(plot = chl_a_plot, filename = here("Output", "chl_a_plot.png"), width = 9, height = 6)

# ANOVA for Chl a content and treatment type #
chla_gen_trtmt_model <- lmer(chla_ug_cm2 ~ TREATMENT + (1|GENOTYPE) +(1|TANK_NUM), data=chl_full_filtered)
check_model(chla_gen_trtmt_model)
summary(chla_gen_trtmt_model)
anova(chla_gen_trtmt_model) # treatment not significant for explaining any variance seen in chla content


chla_norm_summary <- chl_full_filtered %>%
  group_by(TREATMENT) %>%
  summarise(mean_chla_norm = mean(chla_ug_cm2, na.rm = TRUE),
            se_chla_norm = sd(chla_ug_cm2, na.rm = TRUE) / sqrt(n()))

write_csv(chl_full_filtered, here("Data", "Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_chl_full_data.csv"))
