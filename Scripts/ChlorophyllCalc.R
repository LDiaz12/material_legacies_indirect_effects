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
  
## Combine both plates together ## 
plates_comb <- rbind(PlateData1, PlateData2)
plates_comb <- plates_comb[-c(1)]

# metadata
MetaData1<-read_csv(here("Data","Data_Raw","Chl_Content", "Chl_Files", "Metadata1.csv")) 
MetaData2<-read_csv(here("Data","Data_Raw","Chl_Content", "Chl_Files", "Metadata2.csv"))
MetaData2 <- MetaData2 %>% 
  filter(!CORAL_NUM == "EMPTY") ## filter out empty wells 
## Combine plate metadata together ##
metadata_comb <- bind_rows(MetaData1, MetaData2)

##### Analysis #####

# Bring plate and metadata dataframes together #
full_data <- bind_cols(metadata_comb,plates_comb)
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
          Chla_norm = Chla - Chla[CORAL_NUM == "BLANK"], # subtract chl value from blank normalize
          Chlc_norm = Chlc - Chlc[CORAL_NUM == "BLANK"],
          Chl_total = Chla_norm + Chlc_norm)
# after this code, one of the blank's is showing up as having values
# is this because there are two blanks? Should I only normalize to one of them? 

Data_norm2 <- full_join(full_data %>%
                          select(CORAL_NUM, GENOTYPE, TREATMENT, TANK_NUM, CORAL_TANK_NUM), Data_norm, 
                        relationship = "many-to-many")

## Then normalize to SA ##
sa <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA_calculated.csv"))

Data_norm2 <- Data_norm2 %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, TANK_NUM, CORAL_TANK_NUM, Chla_norm, Chlc_norm, Chl_total)

## BELOW ONLY FOR METADATA FILE WRITING ##

sa <- sa %>%
  mutate(CORAL_NUM = as.character(CORAL_NUM))

chl_full <- Data_norm2 %>%
  right_join(sa) %>% 
  select(-c(CORALID, weight1_g, weight2_g, weight_of_wax_g, date))

chl_full <- chl_full %>%                
  mutate(chla.ug.cm2 = Chla_norm * 1 / SA_cm_2, 
         chlc2.ug.cm2 = Chlc_norm * 1 / SA_cm_2) #add in two new columns - chl a normalized to SA and chl c normalized to SA

 
chl_initial <- chl_full %>%
  filter(TREATMENT == "Pre") %>%
  select(c(GENOTYPE, chla.ug.cm2)) %>%
  rename(initial_chla = chla.ug.cm2)

chl_data_full <- chl_full %>%
  left_join(chl_initial)

ggplot(chl_data_full %>%
         filter(!TREATMENT == "Pre")) +
  geom_point(aes(x = initial_chla, y = chla.ug.cm2, color = TREATMENT))

#write_csv(chl_data_full, here("Data", "Data_Raw", "Chl_Content", "Chl_Files", "MO24BEAST_chl_full_data.csv"))

## PLOTS ## 
chl_full_filtered$TREATMENT <- factor(chl_full_filtered$TREATMENT, levels = c("Pre", "Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

chl_total_plot <- ggplot(chl_full_filtered, aes(x=TREATMENT, y=Chl_total, fill=TREATMENT)) + 
  geom_boxplot() + 
  labs(x="Treatment", y="Total Chl (µg/cm2)") + 
  scale_fill_manual(values=c("Control"="blue","Algae_Dom"="darkgreen",
                             "Coral_Dom" = "coral","Rubble_Dom" = "tan", "Pre" = "hotpink"), name="Treatment") +
  theme_bw(base_size=14) +
  geom_point(position = position_dodge(width=0.4))
chl_total_plot # this is just looking at total chl for visualization - not normalized to SA


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
  ggplot(aes(x = TREATMENT, y = chla.ug.cm2, color = TREATMENT)) + # color by treatment for talks
  labs(x = "Treatment", y = expression(bold("Chlorophyll a" ~ (µg ~ cm^-2)))) +
  scale_x_discrete(labels=c("Pre" = "Pre-Experiment", "Algae_Dom" = "Algae-Dominated", "Control" = "Control",
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
  scale_color_manual(values = c("Pre" = "hotpink", "Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                     "Rubble_Dom" = "tan"))
chl_a_plot
#ggsave(plot = chl_a_plot, filename = here("Output", "chl_a_plot.png"), width = 9, height = 6)

# ANOVA for Chl a content and treatment type #
chla_gen_trtmt_model <- lmer(chla.ug.cm2 ~ TREATMENT + (1|GENOTYPE) +(1|TANK_NUM), data=chl_full_filtered)
check_model(chla_gen_trtmt_model)
summary(chla_gen_trtmt_model)
anova(chla_gen_trtmt_model) # treatment not significant for explaining any variance seen in chla content
##############
# chl c plot  # not using Chl c for presentations
chl_full_filtered$TREATMENT <- factor(chl_full_filtered$TREATMENT, levels = c("Pre", "Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

chl_c_plot <- chl_full_filtered %>%
  ggplot(aes(x = TREATMENT, y = chlc2.ug.cm2, color = TREATMENT)) + 
  labs(x = "Treatment", y = "chlorophyll c (µg/cm2)", title = "Chlorophyll c Content by Treatment") +
  geom_jitter(width = 0.1) +
  theme(axis.title = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_line(color = "gray")) +
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1), 
               geom = "errorbar", color = "black", width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", size = 3.5, color = "black") +
  scale_color_manual(values = c("Pre" = "hotpink", "Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
chl_c_plot
#ggsave(plot = chl_c_plot, filename = here("Output", "chl_c_plot.png"), width = 9, height = 6)

# ANOVA for Chl c content and treatment type #
chlc_gen_trtmt_model <- lmer(chlc2.ug.cm2 ~ TREATMENT + (1|GENOTYPE) + (1|TANK_NUM), data=chl_full_filtered)
check_model(chlc_gen_trtmt_model) 
summary(chlc_gen_trtmt_model)
anova(chlc_gen_trtmt_model) # no sig effect of Treatment on explaining variance in chlc2
#################

# Calculate change in chl a and c content from initial starting corals 
chl_initial <- chl_full_filtered %>%
  filter(TREATMENT == "Pre") %>% 
  select(GENOTYPE, CORAL_NUM, Chla_norm_initial = Chla_norm, Chlc_norm_initial = Chlc_norm,
         Chl_total_initial = Chl_total, chla.ug.cm2.initial = chla.ug.cm2, 
         chlc2.ug.cm2.initial = chlc2.ug.cm2)

full_data2 <- chl_full_filtered %>%
  full_join(chl_initial, by = "GENOTYPE", "CORAL_NUM") %>%
  filter(TREATMENT !="Pre") %>%
  mutate(Chla_percent_change = ((Chla_norm - Chla_norm_initial)/Chla_norm_initial) * 100,
         Chlc_percent_change = ((Chlc_norm - Chlc_norm_initial)/Chlc_norm_initial) * 100) %>%
  mutate(TREATMENT = factor(TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))) %>%
  select(-(CORAL_NUM.y))

## Change chl a plot ##
# calculate means of the differences in chl a per treatment
chla_summary <- full_data2 %>%
  group_by(TREATMENT) %>%
  summarise(mean_Chla_percent_change = mean(Chla_percent_change, na.rm = TRUE),
            se_Chla_percent_change = sd(Chla_percent_change, na.rm = TRUE)/sqrt(n()))

chla_norm_summary <- full_data2 %>%
  group_by(TREATMENT) %>%
  summarise(mean_Chla_norm = mean(Chla_norm, na.rm = TRUE),
            se_Chla_norm = sd(Chla_norm, na.rm = TRUE) / sqrt(n()))

## plot change in chl a content from final minus initial corals per treatment with errors 
chla_change_plot <- full_data2 %>%
  ggplot(aes(x = TREATMENT, y = Chla_percent_change, color = TREATMENT)) +
  geom_hline(yintercept = 0) +
  labs(x = "Treatment", y = "Chl a Percent Change") +
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
  #geom_text(data = chla_summary, 
            #aes(x = TREATMENT, y = mean_Chla_percent_change, 
                #label = paste0("", round(mean_Chla_percent_change, 2), "±", round(se_Chla_percent_change, 2))),
           # vjust = -1, hjust = 1.5, color = "black", size = 4)

chla_change_plot
ggsave(plot = chla_change_plot, filename = here("Output", "chla_change_plot.png"), width = 9, height = 6)

delta_chla_model <- lmer(Chla_percent_change ~ TREATMENT + (1|GENOTYPE) +(1|TANK_NUM), data=full_data2)
check_model(delta_chla_model)
summary(delta_chla_model)
anova(delta_chla_model) # no sig effect of treatment on percent change of chla 

#########################
## Change Chl c2 plot ##
chlc_change_plot <- full_data2 %>%
  ggplot(aes(x = TREATMENT, y = Chlc_percent_change, color = TREATMENT)) +
  geom_hline(yintercept = 0) +
  labs(x = "Treatment", y = "Change in Chl c2 Content", title = "Change in Chlorophyll c2 Content by Treatment") +
  geom_jitter(width = 0.1) +
  theme(axis.title = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_line(color = "gray")) +
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    
               geom = "errorbar", color = "black", width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", size = 3.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "#E31A1C", "Control" = "green4", "Coral_Dom" = "dodgerblue2",
                                "Rubble_Dom" = "#6A3D9A"))
chlc_change_plot
ggsave(plot = chlc_change_plot, filename = here("Output", "chlc_change_plot.png"), width = 9, height = 6)
########################


### Combine chlorophyll data and tissue biomass data ###
#chl_full$CORAL_NUM <- as.numeric(chl_full$CORAL_NUM)
afdw_data <- read_csv(here("Data", "Data_Raw", "Growth", "coral_mean_biomass_calculated.csv"))
#afdw_data$CORAL_NUM <- as.numeric(afdw_data$CORAL_NUM)
sa$CORAL_NUM <- as.numeric(sa$CORAL_NUM)
afdw_sa <- right_join(afdw_data, sa)

afdw_sa_noNA <- afdw_sa %>% 
  drop_na()

full_data2$CORAL_NUM.x <- as.numeric(full_data2$CORAL_NUM.x)

chl_biomass_data <- full_data2 %>%
  full_join(afdw_sa_noNA) %>%
  select(TANK_NUM, TREATMENT, GENOTYPE, CORAL_NUM, Chla_norm, Chla_norm_initial, Chla_percent_change,
         mean_tissue_biomass) %>%
  drop_na()

chl_biomass_data$GENOTYPE <- as.character(chl_biomass_data$GENOTYPE)
### regression of chl a and mean tissue biomass ###
chl_biomass_model <- lm(Chla_norm ~ mean_tissue_biomass, data = chl_biomass_data) ## does the amount of chla depend on coral tissue biomass?
check_model(chl_biomass_model)
summary(chl_biomass_model) # small sig effect. chl a and mean tissue biomass slightly related

chl_biomass_reg <- chl_biomass_data %>%
  ggplot(aes(x=mean_tissue_biomass, y=Chla_norm, color = TREATMENT)) +
  geom_point() +
  facet_wrap(~TREATMENT) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_smooth(method = "lm", formula = y~x) + 
  stat_regline_equation(label.x = 1e-04, label.y = 6) + 
  stat_cor(label.x = 1e-04, label.y = 5.5) +
  labs(y = expression(bold("Chla Content" ~ (µg ~ cm^-2))), x= expression(bold("Mean Coral Tissue Biomass" ~ (g ~ mL^-1 ~ cm^-2)))) +
  scale_color_manual (values = c("Algae_Dom" = "darkgreen", "Control" = "blue", 
                                 "Coral_Dom" = "coral", "Rubble_Dom" = "tan"))
chl_biomass_reg
#ggsave(plot = chl_biomass_reg, filename = here("Output", "chl_biomass_reg.png"), width = 9, height = 7)

### Combine chl data and carb chem data ### 
pH_clean <- read_csv(here("Data", "Chemistry", "Cleaned_pH_Data_FULL.csv"))

pH_plotdata<- pH_clean %>%
  group_by(TREATMENT, TANK_NUM) %>%
  summarize(pH_rangemean = mean(pH_range, na.rm = TRUE),
            pH_rangese = sd(pH_range, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH_dailymean, na.rm = TRUE),
            pH_se = sd(pH_dailymean, na.rm = TRUE)/sqrt(n()))

chl_chem_data <- full_data2 %>%
  select(TANKID, TREATMENT, Chla_norm, Chla_diff, GENOTYPE) %>%
  #group_by(TANKID, TREATMENT) %>%
  #summarize(Chla_diff_mean = mean(Chla_diff, na.rm = TRUE),
            #Chla_norm_mean = mean(Chla_norm, na.rm = TRUE)) %>%
  filter(Chla_norm < 6) %>%
  full_join(pH_plotdata)

## plot chla normalized and mean pH in each treatment ##
chla_meanpH_plot <- chl_chem_data %>%
  ggplot(aes(x = pH_mean, y = Chla_norm)) + 
  geom_point(aes(color = TREATMENT)) +
  labs(y = expression(bold("Chla Content" ~ (µg ~ cm^-2))), x= "Mean pH") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_smooth(method = "lm", formula = y~x)
chla_meanpH_plot
ggsave(plot = chla_meanpH_plot, filename = here("Output", "chla_meanpH_plot.png"), width = 9, height = 7)

chl_meanpH_model <- lmer(Chla_norm ~ pH_mean + (1|TANKID), data = chl_chem_data)
anova(chl_pH_model)
summary(chl_pH_model)


## we see a significant relationship between increasing mean pH and increasing chla content 
## for the most part, the control and coral dominated treatments have a mean pH between 8.03 - 8.045
##   and varying chl a content 
## the rubble and almost all of the algae-dominated corals have a mean pH between 8.05 and 8.06
##   the highest chla content in the response corals is only seen in the algae-dominated treatments 


chl_chem_data %>%
  ggplot(aes(x = pH_rangemean, y = Chla_norm)) + 
  geom_point(aes(color = TREATMENT)) +
  geom_smooth(method = "lm")

chl_chem_data %>%
  ggplot(aes(x = pH_mean, y = Chla_diff)) + 
  geom_point(aes(color = TREATMENT)) +
  geom_smooth(method = "lm")

chl_pH_model <- lmer(Chla_norm ~ pH_mean + (1|TANKID), data = chl_chem_data)
anova(chl_pH_model)
summary(chl_pH_model)



chl_chem_plotdata <- chl_chem_data_clean %>%
  group_by(TANKID, TREATMENT, Chla_norm) %>%
  summarize(pH_mean = mean(pH, na.rm = TRUE),
            pH_se = sd(pH, na.rm = TRUE)/sqrt(n()))

chl_pH_model <- lm(Chla_norm ~ pH_mean, chl_chem_plotdata)
plot(chl_pH_model)
summary(chl_pH_model)

chl_pH_reg <- chl_chem_plotdata %>%
  ggplot(aes(x=pH_mean, y=Chla_norm, color = TREATMENT)) +
  geom_point() +
  facet_wrap(~ TREATMENT) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_smooth(method = "lm", formula = y~x) + 
  labs(y = expression(bold("Chla Content" ~ (µg ~ cm^-2))), x="Mean pH")
chl_pH_reg
ggsave(plot = chl_pH_reg, filename = here("Output", "chl_pH_reg.png"), width = 9, height = 7)

###########################
## Change chl c stats ## 
delta_chlc_model_r <- lmer(Chlc_diff~TREATMENT + (1|GENOTYPE), data=full_data2)
plot(delta_chlc_model_r) 
qqp(residuals(delta_chlc_model_r), "norm")
summary(delta_chlc_model_r)
anova(delta_chlc_model_r)
# without random effects
delta_chlc_model <- lm(Chlc_diff ~ TREATMENT, data=full_data2)
qqp(residuals(delta_chlc_model), "norm")
summary(delta_chlc_model)
anova(delta_chlc_model)

anova(delta_chlc_model_r, delta_chlc_model)
