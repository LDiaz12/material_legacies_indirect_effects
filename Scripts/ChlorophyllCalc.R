### Script to read in chlorophyll data from the Synergy 96 well plate in the Molecular Lab ####
### By Nyssa Silbiger ###
### Created on 6/14/2023 #####
### Edited by Laurel Diaz 9/12/24 ###

### load libraries ######
library(tidyverse)
library(lubridate)
library(here)
library(janitor)

### read in data for plate 1 ###
PlateData1<-read_csv(here("Data","Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_Chl_Run1_Plate1.csv"), skip = 39) #skips first 39 lines

# metadata
MetaData<-read_csv(here("Data","Data_Raw","Chl_Content", "Chl_Files", "Metadata1.csv")) 

##### Analysis ####

# Bring together the plate and metadata
Data_combined<-MetaData %>% 
  left_join(PlateData1) %>%
  # Calculate chl from Jeffrey and Humphrey (1975)
  # units in µg/ml
  mutate(Chla = (11.43*(`663`-`750`) - 0.64*(`630` - `750`))*1/0.6,
         Chlc = (27.09*(`630` - `750`) - 3.63*(`663` - `750`)*1/0.6)) 
### Chl a notes ###
# 11.43 is the extinction coefficient of chla at wavelength 663
# 0.64 is the extinction coefficient of chla at wavelength 630

### Chl c2 notes ###
# 27.09 is the extinction coefficient of chl c2 at wavelength 630
# 3.63 is the extinction coefficient of chl c2 at wavelength 663

# 1ml sample
# 0.6 cm path length adjustment

## 1st normalize data to blank ##
Data_norm <- Data_combined %>% 
  group_by(CORAL_NUM)%>%
  ungroup()%>% # ungroup data so I can just look at coral num and chlorophylls
  reframe(CORAL_NUM = CORAL_NUM,
          Chla_norm = Chla - Chla[CORAL_NUM == "BLANK"], # subtract chla value from blank normalize
          Chlc_norm = Chlc - Chlc[CORAL_NUM == "BLANK"],
          Chl_total = Chla_norm + Chlc_norm)
Data_norma <- Data_combined %>% # here I'm rejoining the normalized data to the combined data 
  left_join(Data_norm) %>%
  drop_na() # removing NA values

## Then normalize to SA ##
sa <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA.csv"))
sa$CORAL_NUM <- as.character(sa$CORAL_NUM) # for some reason CORAL_NUM was reading in as double so change to character
chl_1_full <- full_join(Data_norma, sa, relationship = "many-to-many") #full join normalized chl data and surface area data
chl_1_full <- chl_1_full %>%                #many to many relationship bc multiple records associated with multiple tables
  mutate(chla.ug.cm2 = Chla * 1 / SA_cm_2, #need to fix homogenate volume
         chlc2.ug.cm2 = Chlc * 1 / SA_cm_2) %>% #add in two new columns - chl a normalized to SA and chl b normalized to SA
  drop_na() # drop nas because not all SA values are in this first plate

## Plots ## 
chl_total_plot <- ggplot(chl_1_full, aes(x=TREATMENT, y=Chl_total, fill=TREATMENT)) + 
  geom_boxplot() + 
  labs(x="Treatment", y="Total Chl (µg/cm2)") + 
  scale_fill_manual(values=c("Control"="green","Algae_Dom"="red","Coral_Dom" = "blue","Rubble_Dom" = "yellow", "Pre" = "purple"), name="Treatment") +
  theme_bw(base_size=14) +
  geom_point(position = position_dodge(width=0.4))
chl_total_plot # this is just looking at total chl for visualization - NOT normalized
ggsave(plot = chl_total_plot, filename = here("Output", "chl_total_plot.png"), width = 9, height = 6)
## PRE treatment signifies chl measurements from corals BEFORE experiment, so they were not
## put into a treatment

# chl a plot for plate 1 #
chl_a_plot <- chl_1_full %>%
  ggplot(aes(x = TREATMENT, y = chla.ug.cm2, color = TREATMENT)) +
  labs(x = "Treatment", y = "chlorophyll a (µg/cm2)") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun.y = mean, geom = "point", color = "black")   
chl_a_plot # a couple outliers here
ggsave(plot = chl_a_plot, filename = here("Output", "chl_a_plot.png"), width = 9, height = 6)

# chl c plot for plate 1 #
chl_c_plot <- chl_1_full %>%
  ggplot(aes(x = TREATMENT, y = chlc2.ug.cm2, color = TREATMENT)) +
  labs(x = "Treatment", y = "chlorophyll c2 (µg/cm2)") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun.y = mean, geom = "point", color = "black")   
chl_c_plot # what's up with this outlier? 
ggsave(plot = chl_c_plot, filename = here("Output", "chl_c_plot.png"), width = 9, height = 6)

### read in data for plate 2 ###
PlateData2<-read_csv(here("Data","Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_Chl_Run1_Plate2.csv"), skip = 39) #skips first 39 lines

# metadata
MetaData2<-read_csv(here("Data","Data_Raw","Chl_Content", "Chl_Files", "Metadata2.csv")) 

##### Analysis ####

# Bring together the plate and metadata
Data_combined2<-MetaData2 %>% 
  left_join(PlateData2) %>%
  # Calculate chl from Jeffrey and Humphrey (1975)
  # units in µg/ml
  mutate(Chla = (11.43*(`663`-`750`) - 0.64*(`630` - `750`))*1/0.6,
         Chlc = (27.09*(`630` - `750`) - 3.63*(`663` - `750`)*1/0.6)) 
### Chl a notes ###
# 11.43 is the extinction coefficient of chla at wavelength 663
# 0.64 is the extinction coefficient of chla at wavelength 630

### Chl c2 notes ###
# 27.09 is the extinction coefficient of chl c2 at wavelength 630
# 3.63 is the extinction coefficient of chl c2 at wavelength 663

# 1ml sample
# 0.6 cm path length adjustment

## Normalize to blanks ###
Data_combined2 <- Data_combined2 %>% 
  filter(!(CORAL_NUM == "EMPTY"))

Data_norm2 <- Data_combined2 %>% 
  group_by(CORAL_NUM)%>%
  ungroup() %>%
  reframe(CORAL_NUM = CORAL_NUM,
          Chla_norm = Chla - Chla[CORAL_NUM == "BLANK"], 
          Chlc_norm = Chlc - Chlc[CORAL_NUM == "BLANK"],
          Chl_total = Chla_norm + Chlc_norm)
Data_norm2a <- Data_combined2 %>%
  left_join(Data_norm2) %>%
  drop_na()

## Normalize to surface area ##
#sa <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA.csv")) #already read this in
#sa$CORAL_NUM <- as.character(sa$CORAL_NUM) #already did this above
chl_2_full <- full_join(Data_norm2a, sa, relationship = "many-to-many")
chl_2_full <- chl_2_full %>%
  mutate(chla.ug.cm2 = Chla * 1 / SA_cm_2, #need to fix homogenate volume
         chlc2.ug.cm2 = Chlc * 1 / SA_cm_2) %>%
  drop_na()

## Plots ##
# Total chl plot for plate 2 - visualization only  # 
chl_total_plot_2 <- ggplot(chl_2_full, aes(x=TREATMENT, y=Chl_total, fill=TREATMENT)) + 
  geom_boxplot() + 
  labs(x="Treatment", y="Total Chl (µg/cm2)") + 
  scale_fill_manual(values=c("Control"="green","Algae_Dom"="red","Coral_Dom" = "blue", "Pre" = "purple", "Rubble_Dom" = "yellow"), name="Treatment") +
  theme_bw(base_size=14) +
  geom_point(position = position_dodge(width=0.4))
chl_total_plot_2
ggsave(plot = chl_total_plot_2, filename = here("Output", "chl_total_plot_2.png"), width = 9, height = 6)

# chl a plot for plate 2 #
chl_a_plot_2 <- chl_2_full %>%
  ggplot(aes(x = TREATMENT, y = chla.ug.cm2, color = TREATMENT)) +
  labs(x = "Treatment", y = "chlorophyll a (µg/cm2)") +
  geom_jitter(width = 0.1) +                                          
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun.y = mean, geom = "point", color = "black")   
chl_a_plot_2
ggsave(plot = chl_a_plot_2, filename = here("Output", "chl_a_plot_2.png"), width = 9, height = 6)

# chl c plot for plate 2 #
chl_c_plot_2 <- chl_2_full %>%
  ggplot(aes(x = TREATMENT, y = chlc2.ug.cm2, color = TREATMENT)) +
  labs(x = "Treatment", y = "chlorophyll c2 (µg/cm2)") +
  geom_jitter(width = 0.1) +                                           
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),   
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun.y = mean, geom = "point", color = "black")   
chl_c_plot_2 # one weird outlier in rubble dom
ggsave(plot = chl_c_plot_2, filename = here("Output", "chl_c_plot_2.png"), width = 9, height = 6)
