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

### read in plate data ###
PlateData1<-read_csv(here("Data","Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_Chl_Run1_Plate1.csv"), skip = 39) #skips first 39 lines
PlateData2<-read_csv(here("Data","Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_Chl_Run1_Plate2.csv"), skip = 39)
PlateData2 <- PlateData2[-(24:96),] # deleting the rows with empty wells
  
## Combine both plates together ## 
plates_comb <- rbind(PlateData1, PlateData2)
plates_comb <- plates_comb[-c(1)]

# metadata
MetaData1<-read_csv(here("Data","Data_Raw","Chl_Content", "Chl_Files", "Metadata1.csv")) 
MetaData2<-read_csv(here("Data","Data_Raw","Chl_Content", "Chl_Files", "Metadata2.csv"))
MetaData2 <- MetaData2 %>% 
  filter(!CORAL_NUM == "EMPTY") ## filter out empty wells 
## Combine plate metadata together ##
metadata_comb <- rbind(MetaData1, MetaData2)

##### Analysis #####

# Bring plate and metadata data frames together #
full_data <- cbind(metadata_comb,plates_comb)
full_data <- full_data %>% 
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

### we subtract from wavelength of 750 because this is a correction for the turbidity of the sample
# 1ml sample
# 0.6 cm path length adjustment

## 1st normalize data to blank ##
Data_norm <- full_data %>% 
  group_by(CORAL_NUM)%>%
  ungroup()%>% # ungroup data so I can just look at coral num and chlorophylls
  reframe(CORAL_NUM = CORAL_NUM,
          Chla_norm = Chla - Chla[CORAL_NUM == "BLANK"], # subtract chla value from blank normalize
          Chlc_norm = Chlc - Chlc[CORAL_NUM == "BLANK"],
          Chl_total = Chla_norm + Chlc_norm)
Data_norm2 <- cbind(full_data, Data_norm)
Data_norm2 <- Data_norm2[-c(11)]

## Then normalize to SA ##
sa <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA.csv"))
sa$CORAL_NUM <- as.character(sa$CORAL_NUM) # for some reason CORAL_NUM was reading in as double so change to character
chl_full <- full_join(Data_norm2, sa) #full join normalized chl data and surface area data
chl_full <- chl_full %>%                
  mutate(chla.ug.cm2 = Chla * 1 / SA_cm_2, 
         chlc2.ug.cm2 = Chlc * 1 / SA_cm_2) %>% #add in two new columns - chl a normalized to SA and chl b normalized to SA
  drop_na() # drop nas for now but not sure why I'm missing some from the chl?

## PLOTS ## 
chl_total_plot <- ggplot(chl_full, aes(x=TREATMENT, y=Chl_total, fill=TREATMENT)) + 
  geom_boxplot() + 
  labs(x="Treatment", y="Total Chl (µg/cm2)") + 
  scale_fill_manual(values=c("Control"="green","Algae_Dom"="red","Coral_Dom" = "blue","Rubble_Dom" = "yellow", "Pre" = "purple"), name="Treatment") +
  theme_bw(base_size=14) +
  geom_point(position = position_dodge(width=0.4))
chl_total_plot # this is just looking at total chl for visualization - not normalized to SA
ggsave(plot = chl_total_plot, filename = here("Output", "chl_total_plot.png"), width = 9, height = 6)
## PRE treatment signifies chl measurements from corals BEFORE experiment, so they were not put into a treatment

# chl a plot  #
chl_a_plot <- chl_full %>%
  ggplot(aes(x = TREATMENT, y = chla.ug.cm2, color = TREATMENT)) +
  labs(x = "Treatment", y = "chlorophyll a (µg/cm2)") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun.y = mean, geom = "point", color = "black")   
chl_a_plot # a couple outliers here
ggsave(plot = chl_a_plot, filename = here("Output", "chl_a_plot.png"), width = 9, height = 6)

# ANOVA for Chl a content and treatment type #
chla_gen_trtmt_model <- lmer(Chla_norm~TREATMENT + (1|GENOTYPE) +(1|TANKID) + (1|GENOTYPE:TREATMENT) + (1|TANKID:TREATMENT), data=chl_full)
plot(chla_gen_trtmt_model)
# Bonferroni test to identify outliers with undue influence
# Also check with Cook's D
influencePlot(chla_gen_trtmt_model)
# 4 points are greater than 0.1 and are potential outliers
outlierTest(chla_gen_trtmt_model) #Bonferroni p > 0.05
qqp(residuals(chla_trtmt_model), "norm")
summary(chla_gen_trtmt_model)
anova(chla_gen_trtmt_model)
# without genotype as random effect
chla_trtmt_model <- lm(Chla_norm ~ TREATMENT, data=chl_full)
plot(chla_trtmt_model)
qqp(residuals(chla_trtmt_model), "norm")
summary(chla_trtmt_model)
anova(chla_trtmt_model)

# chl c plot  #
chl_c_plot <- chl_full %>%
  ggplot(aes(x = TREATMENT, y = chlc2.ug.cm2, color = TREATMENT)) +
  labs(x = "Treatment", y = "chlorophyll c2 (µg/cm2)") +
  geom_jitter(width = 0.1) +                                            
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun.y = mean, geom = "point", color = "black")   
chl_c_plot # what's up with this outlier? 
ggsave(plot = chl_c_plot, filename = here("Output", "chl_c_plot.png"), width = 9, height = 6)

# ANOVA for Chl c content and treatment type #
chlc_gen_trtmt_model <- lmer(Chlc_norm~TREATMENT + (1|GENOTYPE) + (1|TANKID) + (1|GENOTYPE:TREATMENT) + (1|TANKID:TREATMENT), data=chl_full)
plot(chlc_gen_trtmt_model) #yikes
influencePlot(chlc_gen_trtmt_model) # 4 potential outliers with undue influence
outlierTest(chlc_gen_trtmt_model) #Bonferroni p < 0.05
qqp(residuals(chlc_gen_trtmt_model), "norm")
summary(chlc_gen_trtmt_model)
anova(chlc_gen_trtmt_model)
# without genotype as random effect
chlc_trtmt_model <- lm(Chlc_norm ~ TREATMENT, data=chl_full)
qqp(residuals(chlc_trtmt_model), "norm")
summary(chlc_trtmt_model)
anova(chlc_trtmt_model)
