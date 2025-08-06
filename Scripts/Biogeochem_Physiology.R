library(tidyverse)
library(here)
library(ggridges)
library(agricolae)
library(lme4)
library(lmerTest)
library(moments)
library(performance)
library(ggpubr)
library(emmeans)
library(cowplot)
#install.packages("sjPlot")
library(sjPlot)
library(sjmisc)
#install.packages("effects")
library(effects)
library(sjstats)
library(patchwork)


## ANALYSIS OF BIOGEOCHEM IMPACTING CORAL PHYSIOLOGY ## 

# load in metadata sheet # 
metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))
raw_chem <- read_csv(here("Data", "Chemistry", "Raw_Chem_Data.csv"))

metadata_raw_chem <- metadata %>%
  left_join(raw_chem)

metadata_raw_chem_means <- metadata_raw_chem %>%
  group_by(TANK_NUM) %>%
  mutate(mean_endos = mean(endo_per_cm2, na.rm = TRUE),
         mean_chl = mean(chla_ug_cm2, na.rm = TRUE), 
         mean_biomass = mean(mean_tissue_biomass, na.rm = TRUE), 
         mean_GP = mean(GP, na.rm = TRUE),
         mean_R = mean(R, na.rm = TRUE))

####### ANOVAS OF PHYSIO PARAMS PER TREATMENT #######
metadata_raw_chem_means$TREATMENT <- factor(metadata_raw_chem_means$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

# ENDOS #
endos_plot <- metadata1 %>% 
  ggplot(aes(x = TREATMENT, y = endo_per_cm2, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Endosymbiont Density" ~ (cells ~ x10^6 ~ cm^-2)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata, aes(x = TREATMENT, y = endo_per_cm2), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
endos_plot
#ggsave(plot = endos_plot, filename = here("Output", "Biogeochem_Physio", "endos_plot.png"), width = 14, height = 10)

endos_model <- lmer(endo_per_cm2 ~ TREATMENT + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data= metadata1)
check_model(endos_model)
summary(endos_model)
anova(endos_model) # non significant. p = 0.65

# CHL A #
chla_plot <- metadata1 %>% 
  ggplot(aes(x = TREATMENT, y = chla_ug_cm2, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Chlorophyll-a Content" ~ (µg ~ cm^-2)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata1, aes(x = TREATMENT, y = chla_ug_cm2), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
chla_plot
#ggsave(plot = chla_plot, filename = here("Output", "Biogeochem_Physio", "chla_plot.png"), width = 14, height = 10)

chla_model <- lmer(chla_ug_cm2 ~ TREATMENT + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata1)
check_model(chla_model)
summary(chla_model)
anova(chla_model) 


# MEAN TISSUE BIOMASS # 
tissuebiomass_plot <- metadata1 %>% 
  ggplot(aes(x = TREATMENT, y = mean_tissue_biomass, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Mean Tissue Biomass" ~ (mg ~ cm^-2)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata1, aes(x = TREATMENT, y = mean_tissue_biomass), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
tissuebiomass_plot
#ggsave(plot = tissuebiomass_plot, filename = here("Output", "Biogeochem_Physio", "tissuebiomass_plot.png"), width = 14, height = 10)

tissuebiomass_model <- lmer(mean_tissue_biomass ~ TREATMENT + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata1)
check_model(tissuebiomass_model)
summary(tissuebiomass_model)
anova(tissuebiomass_model) 


# R #
R_plot <- metadata1 %>% 
  ggplot(aes(x = TREATMENT, y = R, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Respiration Rate" ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata1, aes(x = TREATMENT, y = R), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
R_plot
#ggsave(plot = R_plot, filename = here("Output", "Biogeochem_Physio", "R_plot.png"), width = 14, height = 10)

R_model <- lmer(R ~ TREATMENT + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data=metadata1)
check_model(R_model)
summary(R_model)
anova(R_model) # significant. p < 0.01. p = 0.0015


# GP # 
GP_plot <- metadata1 %>% 
  ggplot(aes(x = TREATMENT, y = GP, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Gross Photosynthesis" ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata1, aes(x = TREATMENT, y = GP), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
GP_plot
#ggsave(plot = GP_plot, filename = here("Output", "Biogeochem_Physio", "GP_plot.png"), width = 10, height = 10)

GP_model <- lmer(GP ~ TREATMENT + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data=metadata1)
check_model(GP_model)
summary(GP_model)
anova(GP_model) # significant. p < 0.05. p = 0.013

physio_params_patch <- (endos_plot + chla_plot + tissuebiomass_plot)/(R_plot + GP_plot) + plot_annotation(tag_levels = "a")
physio_params_patch
ggsave(plot = physio_params_patch, filename = here("Output", "Biogeochem_Physio", "physio_params_patch.png"), width = 11, height = 11)

######## BIOGEOCHEM IMPACTS ON CHL-A ########

# regression of chl-a and mean pH # 
chla_meanpH_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_pH, y = mean_chl)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x, color = "black") + 
  labs(x = "pH", y = "Chlorophyll-a ug/cm2") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
chla_meanpH_plot
#ggsave(plot = chla_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "chla_meanpH.png"), width = 14, height = 10)

chla_meanpH_model <- lm(mean_chl ~ grand_mean_pH, data = metadata_raw_chem_means)
check_model(chla_meanpH_model)
summary(chla_meanpH_model) #p < 0.001

# regression of chl-a and mean DOC # 
chla_meanDOC_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_DOC, y = mean_chl)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x, color = "black") + 
  labs(x = "DOC (NPOC uM)", y = "Chlorophyll-a ug/cm2") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  #stat_regline_equation(label.x = 160, label.y = 0.14, size = 8) + 
  #stat_cor(label.x = 160, label.y = 0.13, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
chla_meanDOC_plot
#ggsave(plot = chla_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "chla_meanDOC.png"), width = 14, height = 10)

chla_meanDOC_model <- lm(mean_chl ~ grand_mean_DOC, data = metadata_raw_chem_means)
check_model(chla_meanDOC_model)
summary(chla_meanDOC_model) #sig. p = 0.004 

######## BIOGEOCHEM IMPACTS ON ENDO DENSITY #######

# regression of endos and mean pH #
endo_meanpH_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_pH, y = mean_endos)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean pH", y = "Endosymbiont Density (counts/cm2)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  #stat_regline_equation(label.x = 8.043, label.y = 1.3, size = 8) + 
  #stat_cor(label.x = 8.043, label.y = 1.2, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
endo_meanpH_plot
#ggsave(plot = endo_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "endo_meanpH.png"), width = 14, height = 10)

endo_meanpH_model <- lm(mean_endos ~ grand_mean_pH, data = metadata_raw_chem_means)
check_model(endo_meanpH_model)
summary(endo_meanpH_model)

# regression of endos and mean DOC #
endo_meanDOC_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_DOC, y = mean_endos)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean DOC", y = "Endosymbiont Density (counts/cm2)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  #stat_regline_equation(label.x = 150, label.y = 1.3, size = 8) + 
  #stat_cor(label.x = 150, label.y = 1.2, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
endo_meanDOC_plot
#ggsave(plot = endo_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "endo_meanDOC.png"), width = 14, height = 10)

endo_meanDOC_model <- lm(mean_endos ~ grand_mean_DOC, data = metadata_raw_chem_means)
check_model(endo_meanDOC_model)
summary(endo_meanDOC_model)

####### BIOGEOCHEM IMPACTS ON TISSUE BIOMASS ### #####

# regression of tissue biomass and mean pH # 
biomass_meanpH_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_pH, y = mean_biomass)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean pH", y = "Mean Tissue Biomass") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  #stat_regline_equation(label.x = 8.04, label.y = 4e-04, size = 8) + 
  #stat_cor(label.x = 8.04, label.y = 3.5e-04, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
biomass_meanpH_plot
#ggsave(plot = biomass_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "biomass_meanpH.png"), width = 14, height = 10)

biomass_meanpH_model <- lm(mean_biomass ~ grand_mean_pH, data = metadata_raw_chem_means)
check_model(biomass_meanpH_model)
summary(biomass_meanpH_model)

# regression of tissue biomass and mean DOC # 
biomass_meanDOC_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_DOC, y = mean_biomass)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean DOC", y = "Mean Tissue Biomass") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  #stat_regline_equation(label.x = 150, label.y = 4e-04, size = 8) + 
  #stat_cor(label.x = 150, label.y = 3.5e-04, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
biomass_meanDOC_plot
#ggsave(plot = biomass_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "biomass_meanDOC.png"), width = 14, height = 10)

biomass_meanDOC_model <- lm(mean_biomass ~ grand_mean_DOC, data = metadata_raw_chem_means)
check_model(biomass_meanDOC_model)
summary(biomass_meanDOC_model)

##### BIOGEOCHEM IMPACTS ON RESPIRATION RATE #####

# R and mean pH #
R_meanpH_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_pH, y = mean_R)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean pH", y = "Respiration Rate") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  #stat_regline_equation(label.x = 8.03, label.y = 0.07, size = 8) + 
  #stat_cor(label.x = 8.03, label.y = 0.065, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
R_meanpH_plot
#ggsave(plot = R_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "R_meanpH.png"), width = 14, height = 10)

R_meanpH_model <- lm(mean_R ~ grand_mean_pH, data = metadata_raw_chem_means)
check_model(R_meanpH_model)
summary(R_meanpH_model)

# R and mean DOC # 
R_meanDOC_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_DOC, y = mean_R)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean DOC", y = "Respiration Rate") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  #stat_regline_equation(label.x = 150, label.y = 0.07, size = 8) + 
  #stat_cor(label.x = 150, label.y = 0.065, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
R_meanDOC_plot
#ggsave(plot = R_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "R_meanDOC.png"), width = 14, height = 10)

R_meanDOC_model <- lm(mean_R ~ grand_mean_DOC, data = metadata_raw_chem_means)
check_model(R_meanDOC_model)
summary(R_meanDOC_model)

##### BIOGEOCHEM IMPACTS ON NET PHOTOSYNTHESIS ##### 

# NP and mean pH # 
NP_meanpH_plot <- metadata %>%
  ggplot(aes(x = pH_mean, y = NP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean pH", y = "Net Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 8.05, label.y = 0.14, size = 5) + 
  stat_cor(label.x = 8.05, label.y = 0.13, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
NP_meanpH_plot
#ggsave(plot = NP_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "NP_meanpH.png"), width = 14, height = 12)

# NP and range pH # 
NP_rangepH_plot <- metadata %>%
  ggplot(aes(x = pH_rangemean, y = NP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in pH", y = "Net Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 0.00, label.y = 0.14, size = 5) + 
  stat_cor(label.x = 0.00, label.y = 0.13, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
NP_rangepH_plot
#ggsave(plot = NP_rangepH_plot, filename = here("Output", "Biogeochem_Physio", "NP_rangepH.png"), width = 14, height = 12)

# NP and mean DOC # 
NP_meanDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_mean, y = NP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean DOC", y = "Net Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 130, label.y = 0.14, size = 5) + 
  stat_cor(label.x = 130, label.y = 0.13, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
NP_meanDOC_plot
#ggsave(plot = NP_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "NP_meanDOC.png"), width = 14, height = 12)

# NP and range DOC # 
NP_rangeDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_rangemean, y = NP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in DOC", y = "Net Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 0, label.y = 0.14, size = 5) + 
  stat_cor(label.x = 0, label.y = 0.13, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
NP_rangeDOC_plot
#ggsave(plot = NP_rangeDOC_plot, filename = here("Output", "Biogeochem_Physio", "NP_rangeDOC.png"), width = 14, height = 12)

##### BIOGEOCHEM IMPACTS ON GROSS PHOTOSYNTHESIS #####

# GP and mean pH # 
GP_meanpH_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_pH, y = mean_GP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean pH", y = "Gross Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  #stat_regline_equation(label.x = 8.045, label.y = 0.20, size = 5) + 
  #stat_cor(label.x = 8.045, label.y = 0.19, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
GP_meanpH_plot
#ggsave(plot = GP_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "GP_meanpH.png"), width = 14, height = 12)

GP_meanpH_model <- lm(mean_GP ~ grand_mean_pH, data = metadata_raw_chem_means)
check_model(GP_meanpH_model)
summary(GP_meanpH_model)

# GP and mean DOC # 
GP_meanDOC_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_DOC, y = mean_GP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean DOC", y = "Gross Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  #stat_regline_equation(label.x = 150, label.y = 0.20, size = 5) + 
  #stat_cor(label.x = 150, label.y = 0.19, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
GP_meanDOC_plot
#ggsave(plot = GP_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "GP_meanDOC.png"), width = 14, height = 12)

GP_meanDOC_model <- lm(mean_GP ~ grand_mean_DOC, data = metadata_raw_chem_means)
check_model(GP_meanDOC_model)
summary(GP_meanDOC_model)


##### CALCULATE AVERAGE CORAL PHYSIO PARAMETERS PER TREATMENT #####

coral_physio_params_treatment <- metadata %>%
  group_by(TREATMENT) %>%
  summarize(avg_endos = mean(endo_per_cm2, na.rm = TRUE),
            endo_error = sd(endo_per_cm2, na.rm = TRUE)/sqrt(n()),
            avg_chl = mean(chla_ug_cm2, na.rm = TRUE),
            chl_error = sd(chla_ug_cm2, na.rm = TRUE)/sqrt(n()),
            avg_biomass = mean(mean_tissue_biomass, na.rm = TRUE), 
            biomass_error = sd(mean_tissue_biomass, na.rm = TRUE)/sqrt(n()),
            avg_R = mean(R, na.rm = TRUE), 
            R_error = sd(R, na.rm = TRUE)/sqrt(n()),
            avg_NP = mean(NP, na.rm = TRUE), 
            NP_error = sd(NP, na.rm = TRUE)/sqrt(n()),
            avg_GP = mean(GP, na.rm = TRUE),
            GP_error = sd(GP, na.rm = TRUE)/sqrt(n()),
            avg_NEC = mean(NEC_mean, na.rm = TRUE), 
            NEC_error = sd(NEC_mean, na.rm = TRUE)/sqrt(n()),
            avg_NEP = mean(NEP_mean, na.rm = TRUE),
            NEP_error = sd(NEP_mean, na.rm = TRUE)/sqrt(n()))

#write_csv(coral_physio_params_treatment, here("Data", "Summary_Files", "coral_physio_params_treatment.csv"))

metadata$TREATMENT <- factor(metadata$TREATMENT, levels = c("Control", "Algae_Dom", 
                                                            "Coral_Dom", "Rubble_Dom"))

endos_treatment <- lmer(endo_per_cm2 ~ TREATMENT + (1|GENOTYPE), data = metadata)
summary(endos_treatment)
anova(endos_treatment)

coral_physio_params_treatment$TREATMENT <- factor(coral_physio_params_treatment$TREATMENT, levels = c("Control", "Algae_Dom",
                                                                                                      "Coral_Dom", "Rubble_Dom"))
avg_endos_plot <- coral_physio_params_treatment %>%
  ggplot(aes(x = TREATMENT, y = avg_endos, color = TREATMENT)) + 
  geom_point(data = metadata, aes(x = TREATMENT, y = endo_per_cm2), alpha = 0.25) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1, color = "black") +
  stat_summary(fun.y = mean, geom = "point", size = 3, color = "black") +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_endos_plot

avg_chl_plot <- coral_physio_params_treatment %>%
  ggplot(aes(x = TREATMENT, y = avg_chl, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_chl_plot

avg_biomass_plot <- coral_physio_params_treatment %>%
  ggplot(aes(x = TREATMENT, y = avg_biomass, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_biomass_plot

avg_R_plot <- coral_physio_params_treatment %>%
  ggplot(aes(x = TREATMENT, y = avg_R, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_R_plot

avg_NP_plot <- coral_physio_params_treatment %>%
  ggplot(aes(x = TREATMENT, y = avg_NP, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_NP_plot

avg_GP_plot <- coral_physio_params_treatment %>%
  ggplot(aes(x = TREATMENT, y = avg_GP, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_GP_plot

##### CALCULATE AVERAGE CORAL PHYSIO PARAMS PER TANK ##### 
coral_physio_params_tank <- metadata %>%
  group_by(TANK_NUM, TREATMENT) %>%
  summarize(avg_endos = mean(endo_per_cm2, na.rm = TRUE),
            avg_chl = mean(chla.ug.cm2, na.rm = TRUE), 
            avg_biomass = mean(mean_tissue_biomass, na.rm = TRUE), 
            avg_R = mean(R, na.rm = TRUE), 
            avg_NP = mean(NP, na.rm = TRUE), 
            avg_GP = mean(GP, na.rm = TRUE),
            avg_NEC = mean(NEC_mean, na.rm = TRUE), 
            avg_NEP = mean(NEP_mean, na.rm = TRUE))
coral_physio_params_tank
#write_csv(coral_physio_params_tank, here("Data", "Summary_Files", "coral_physio_params_tank.csv"))

coral_physio_params_tank$TREATMENT <- factor(coral_physio_params_tank$TREATMENT, levels = c("Control", "Algae_Dom",
                                                                                                      "Coral_Dom", "Rubble_Dom"))
avg_endos_plot <- coral_physio_params_tank %>%
  ggplot(aes(x = TREATMENT, y = avg_endos, color = TREATMENT)) + 
  geom_boxplot() + 
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_endos_plot

avg_chl_plot <- coral_physio_params_tank %>%
  ggplot(aes(x = TREATMENT, y = avg_chl, color = TREATMENT)) + 
  geom_boxplot() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_chl_plot

avg_biomass_plot <- coral_physio_params_tank %>%
  ggplot(aes(x = TREATMENT, y = avg_biomass, color = TREATMENT)) + 
  geom_boxplot() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_biomass_plot

avg_R_plot <- coral_physio_params_tank %>%
  ggplot(aes(x = TREATMENT, y = avg_R, color = TREATMENT)) + 
  geom_boxplot() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_R_plot

avg_NP_plot <- coral_physio_params_tank %>%
  ggplot(aes(x = TREATMENT, y = avg_NP, color = TREATMENT)) + 
  geom_boxplot() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_NP_plot

avg_GP_plot <- coral_physio_params_tank %>%
  ggplot(aes(x = TREATMENT, y = avg_GP, color = TREATMENT)) + 
  geom_boxplot() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_GP_plot
