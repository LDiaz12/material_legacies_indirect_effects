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
metadata <- metadata %>%
  mutate(mean_tissue_biomass_ng = (mean_tissue_biomass)*1000000)

metadata <- metadata %>%
  mutate(chla.ug.cm2 = ifelse(chla.ug.cm2 < 0, NA, chla.ug.cm2),
         chla.ug.cm2 = ifelse(chla.ug.cm2 > 0.09, NA, chla.ug.cm2),
         chla.ug.cm2 = ifelse(TREATMENT == "Coral_Dom" & chla.ug.cm2 > 0.06, NA, chla.ug.cm2),
         chla.ug.cm2 = ifelse(TREATMENT == "Control" & chla.ug.cm2 > 0.06, NA, chla.ug.cm2),
         chla.ug.cm2 = ifelse(TREATMENT == "Rubble_Dom" & chla.ug.cm2 > 0.07, NA, chla.ug.cm2),
         endo_per_cm2 = ifelse(endo_per_cm2 > 1, NA, endo_per_cm2),
         mean_tissue_biomass_ng = ifelse(mean_tissue_biomass_ng > 750, NA, mean_tissue_biomass_ng),
         mean_tissue_biomass_ng = ifelse(TREATMENT == "Coral_Dom" & mean_tissue_biomass_ng > 300, NA, mean_tissue_biomass_ng),
         R = ifelse(R > 0.06, NA, R),
         R = ifelse(TREATMENT == "Control" & R > 0.045, NA, R),
         R = ifelse(TREATMENT == "Control" & R < 0.02, NA, R), 
         R = ifelse(TREATMENT == "Algae_Dom" & R < 0.02, NA, R), 
         R = ifelse(TREATMENT == "Coral_Dom" & R < 0.02, NA, R),
         R = ifelse(TREATMENT == "Coral_Dom" & R > 0.045, NA, R),
         R = ifelse(TREATMENT == "Rubble_Dom" & R < 0.02, NA, R), 
         GP = ifelse(TREATMENT == "Control" & GP > 0.15, NA, GP), 
         GP = ifelse(TREATMENT == "Algae_Dom" & GP > 0.20, NA, GP),
         GP = ifelse(TREATMENT == "Algae_Dom" & GP < 0.05, NA, GP),
         GP = ifelse(TREATMENT == "Coral_Dom" & GP > 0.15, NA, GP),
         GP = ifelse(TREATMENT == "Rubble_Dom" & GP > 0.15, NA, GP),
         GP = ifelse(TREATMENT == "Rubble_Dom" & GP < 0.05, NA, GP))

#NEC = ifelse(TREATMENT == "Rubble_Dom" & NEC > 2, NA, NEC)
####### ANOVAS OF PHYSIO PARAMS PER TREATMENT #######
metadata$TREATMENT <- factor(metadata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

# ENDOS #
endos_plot <- metadata %>% 
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

endos_model <- lmer(endo_per_cm2 ~ TREATMENT + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data=metadata)
check_model(endos_model)
summary(endos_model)
anova(endos_model) # non significant. p = 0.7086

# CHL A #
chla_plot <- metadata %>% 
  ggplot(aes(x = TREATMENT, y = chla.ug.cm2, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Chlorophyll-a Content" ~ (µg ~ cm^-2)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata, aes(x = TREATMENT, y = chla.ug.cm2), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
chla_plot
#ggsave(plot = chla_plot, filename = here("Output", "Biogeochem_Physio", "chla_plot.png"), width = 14, height = 10)

chla_model <- lmer(chla.ug.cm2 ~ TREATMENT + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(chla_model)
summary(chla_model)
anova(chla_model) # p = 0.545
 
emmeans(chla_model, pairwise ~ "TREATMENT", adjust = "Tukey") 
chla_model_noRandom <- lm(chla.ug.cm2 ~ TREATMENT, data = metadata)
HSD.test(chla_model_noRandom, "TREATMENT", console=TRUE)


# MEAN TISSUE BIOMASS # 
tissuebiomass_plot <- metadata %>% 
  ggplot(aes(x = TREATMENT, y = mean_tissue_biomass_ng, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Mean Tissue Biomass" ~ (ng ~ cm^-2)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata, aes(x = TREATMENT, y = mean_tissue_biomass_ng), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
tissuebiomass_plot
#ggsave(plot = tissuebiomass_plot, filename = here("Output", "Biogeochem_Physio", "tissuebiomass_plot.png"), width = 14, height = 10)

tissuebiomass_model <- lmer(mean_tissue_biomass_ng ~ TREATMENT + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(tissuebiomass_model)
summary(tissuebiomass_model)
anova(tissuebiomass_model) # significant. (p < 0.001) p = 0.0009667

emmeans(tissuebiomass_model, pairwise ~ "TREATMENT", adjust = "Tukey") # control and algae sig diff p = 0.03
tissuebiomass_model_noRandom <- lm(mean_tissue_biomass_ng ~ TREATMENT, data = metadata)
HSD.test(tissuebiomass_model_noRandom, "TREATMENT", console=TRUE)


# R #
R_plot <- metadata %>% 
  ggplot(aes(x = TREATMENT, y = R, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Respiration Rate" ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata, aes(x = TREATMENT, y = R), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
R_plot
#ggsave(plot = R_plot, filename = here("Output", "Biogeochem_Physio", "R_plot.png"), width = 14, height = 10)

R_model <- lmer(R ~ TREATMENT + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data=metadata)
check_model(R_model)
summary(R_model)
anova(R_model) # significant. p < 0.01. p = 0.0015

emmeans(R_model, pairwise ~ "TREATMENT", adjust = "Tukey") # control v algae and algae v coral sig
R_model_noRandom <- lm(R ~ TREATMENT, data = metadata)
HSD.test(R_model_noRandom, "TREATMENT", console=TRUE)
# rubble = a, algae = ab, coral = bc, control = c

# GP # 
GP_plot <- metadata %>% 
  ggplot(aes(x = TREATMENT, y = GP, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Gross Photosynthesis" ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata, aes(x = TREATMENT, y = GP), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
GP_plot
#ggsave(plot = GP_plot, filename = here("Output", "Biogeochem_Physio", "GP_plot.png"), width = 10, height = 10)

GP_model <- lmer(GP ~ TREATMENT + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data=metadata)
check_model(GP_model)
summary(GP_model)
anova(GP_model) # significant. p < 0.05. p = 0.013

emmeans(GP_model, pairwise ~ "TREATMENT", adjust = "Tukey") # control and algae sig p = 0.05
GP_model_noRandom <- lm(GP ~ TREATMENT, data = metadata)
HSD.test(GP_model_noRandom, "TREATMENT", console=TRUE)
# algae = a, rubble = a, coral = ab, control = b

physio_params_patch <- (endos_plot + chla_plot + tissuebiomass_plot)/(R_plot + GP_plot) + plot_annotation(tag_levels = "a")
physio_params_patch
#ggsave(plot = physio_params_patch, filename = here("Output", "Biogeochem_Physio", "physio_params_patch.png"), width = 11, height = 11)

######## BIOGEOCHEM IMPACTS ON CHL-A ########

# regression of chl-a and mean pH # 
chla_meanpH_plot <- metadata %>%
  ggplot(aes(x = pH_mean, y = chla.ug.cm2)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  stat_regline_equation(label.x = 8.035, label.y = 0.10, size = 8) + 
  stat_cor(label.x = 8.035, label.y = 0.09, size = 8) +
  labs(x = "Daily Mean pH", y = "Chlorophyll-a ug/cm2") +
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

chla_meanpH_model <- lmer(chla.ug.cm2 ~ pH_mean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(chla_meanpH_model)
summary(chla_meanpH_model)

# regression of chl-a and range in pH # 
chla_rangepH_plot <- metadata %>%
  ggplot(aes(x = pH_rangemean, y = chla.ug.cm2)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  stat_regline_equation(label.x = 0.05, label.y = 0.14, size = 8) + 
  stat_cor(label.x = 0.05, label.y = 0.13, size = 8) +
  labs(x = "Daily Mean Range in pH", y = "Chlorophyll-a ug/cm2") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
chla_rangepH_plot
#ggsave(plot = chla_rangepH_plot, filename = here("Output", "Biogeochem_Physio", "chla_rangepH.png"), width = 14, height = 10)

chla_rangepH_model <- lmer(chla.ug.cm2 ~ pH_rangemean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(chla_rangepH_model)
summary(chla_rangepH_model) # significant at p = 0.04

# regression of chl-a and mean DOC # 
chla_meanDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_mean, y = chla.ug.cm2)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean DOC (NPOC uM)", y = "Chlorophyll-a ug/cm2") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 160, label.y = 0.14, size = 8) + 
  stat_cor(label.x = 160, label.y = 0.13, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
chla_meanDOC_plot
#ggsave(plot = chla_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "chla_meanDOC.png"), width = 14, height = 10)

chla_meanDOC_model <- lmer(chla.ug.cm2 ~ DOC_mean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(chla_meanDOC_model)
summary(chla_meanDOC_model)

# regression of chl-a and mean range of DOC # 
chla_rangeDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_rangemean, y = chla.ug.cm2)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in DOC (NPOC uM)", y = "Chlorophyll-a ug/cm2") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 40, label.y = 0.10, size = 8) + 
  stat_cor(label.x = 40, label.y = 0.09, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
chla_rangeDOC_plot
#ggsave(plot = chla_rangeDOC_plot, filename = here("Output", "Biogeochem_Physio", "chla_rangeDOC.png"), width = 14, height = 10)

chla_rangeDOC_model <- lmer(chla.ug.cm2 ~ DOC_rangemean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(chla_rangeDOC_model)
summary(chla_rangeDOC_model)

######## BIOGEOCHEM IMPACTS ON ENDO DENSITY #######

# regression of endos and mean pH #
endo_meanpH_plot <- metadata %>%
  ggplot(aes(x = pH_mean, y = endo_per_cm2)) + 
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
  stat_regline_equation(label.x = 8.043, label.y = 1.3, size = 8) + 
  stat_cor(label.x = 8.043, label.y = 1.2, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
endo_meanpH_plot
#ggsave(plot = endo_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "endo_meanpH.png"), width = 14, height = 10)

endo_meanpH_model <- lmer(endo_per_cm2 ~ pH_mean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(endo_meanpH_model)
summary(endo_meanpH_model)

#regression of endos and mean range in pH # 
endo_rangepH_plot <- metadata %>%
  ggplot(aes(x = pH_rangemean, y = endo_per_cm2)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in pH", y = "Endosymbiont Density (counts/cm2)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 0.03, label.y = 1.3, size = 8) + 
  stat_cor(label.x = 0.03, label.y = 1.2, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
endo_rangepH_plot
#ggsave(plot = endo_rangepH_plot, filename = here("Output", "Biogeochem_Physio", "endo_rangepH.png"), width = 14, height = 10)

endo_rangepH_model <- lmer(endo_per_cm2 ~ pH_rangemean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(endo_rangepH_model)
summary(endo_rangepH_model)

# regression of endos and mean DOC #
endo_meanDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_mean, y = endo_per_cm2)) + 
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
  stat_regline_equation(label.x = 150, label.y = 1.3, size = 8) + 
  stat_cor(label.x = 150, label.y = 1.2, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
endo_meanDOC_plot
#ggsave(plot = endo_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "endo_meanDOC.png"), width = 14, height = 10)

endo_meanDOC_model <- lmer(endo_per_cm2 ~ DOC_mean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(endo_meanDOC_model)
summary(endo_meanDOC_model)

# regression of endos and daily mean range in DOC #
endo_rangeDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_rangemean, y = endo_per_cm2)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in DOC", y = "Endosymbiont Density (counts/cm2)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 40, label.y = 1.3, size = 8) + 
  stat_cor(label.x = 40, label.y = 1.2, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
endo_rangeDOC_plot
#ggsave(plot = endo_rangeDOC_plot, filename = here("Output", "Biogeochem_Physio", "endo_rangeDOC.png"), width = 14, height = 10)

endo_rangeDOC_model <- lmer(endo_per_cm2 ~ DOC_rangemean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(endo_rangeDOC_model)
summary(endo_rangeDOC_model)


####### BIOGEOCHEM IMPACTS ON TISSUE BIOMASS ### #####

# regression of tissue biomass and mean pH # 
biomass_meanpH_plot <- metadata %>%
  ggplot(aes(x = pH_mean, y = mean_tissue_biomass)) + 
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
  stat_regline_equation(label.x = 8.04, label.y = 4e-04, size = 8) + 
  stat_cor(label.x = 8.04, label.y = 3.5e-04, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
biomass_meanpH_plot
#ggsave(plot = biomass_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "biomass_meanpH.png"), width = 14, height = 10)

biomass_meanpH_model <- lmer(mean_tissue_biomass ~ pH_mean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(biomass_meanpH_model)
summary(biomass_meanpH_model)

# regression of tissue biomass and range in pH # 
biomass_rangepH_plot <- metadata %>%
  ggplot(aes(x = pH_rangemean, y = mean_tissue_biomass)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in pH", y = "Mean Tissue Biomass") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 0.00, label.y = 4e-04, size = 8) + 
  stat_cor(label.x = 0.00, label.y = 3.5e-04, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
biomass_rangepH_plot
#ggsave(plot = biomass_rangepH_plot, filename = here("Output", "Biogeochem_Physio", "biomass_rangepH.png"), width = 14, height = 10)

biomass_rangepH_model <- lmer(mean_tissue_biomass ~ pH_rangemean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(biomass_rangepH_model)
summary(biomass_rangepH_model)

# regression of tissue biomass and mean DOC # 
biomass_meanDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_mean, y = mean_tissue_biomass)) + 
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
  stat_regline_equation(label.x = 150, label.y = 4e-04, size = 8) + 
  stat_cor(label.x = 150, label.y = 3.5e-04, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
biomass_meanDOC_plot
#ggsave(plot = biomass_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "biomass_meanDOC.png"), width = 14, height = 10)

biomass_meanDOC_model <- lmer(mean_tissue_biomass ~ DOC_mean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(biomass_meanDOC_model)
summary(biomass_meanDOC_model)

# regression of tissue biomass and range in DOC # 
biomass_rangeDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_rangemean, y = mean_tissue_biomass)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in DOC", y = "Mean Tissue Biomass") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 40, label.y = 4e-04, size = 8) + 
  stat_cor(label.x = 40, label.y = 3.5e-04, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
biomass_rangeDOC_plot
#ggsave(plot = biomass_rangeDOC_plot, filename = here("Output", "Biogeochem_Physio", "biomass_rangeDOC.png"), width = 14, height = 10)

biomass_rangeDOC_model <- lmer(mean_tissue_biomass ~ DOC_rangemean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(biomass_rangeDOC_model)
summary(biomass_rangeDOC_model)

##### BIOGEOCHEM IMPACTS ON RESPIRATION RATE #####

# R and mean pH #
R_meanpH_plot <- metadata %>%
  ggplot(aes(x = pH_mean, y = R)) + 
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
  stat_regline_equation(label.x = 8.03, label.y = 0.07, size = 8) + 
  stat_cor(label.x = 8.03, label.y = 0.065, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
R_meanpH_plot
#ggsave(plot = R_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "R_meanpH.png"), width = 14, height = 10)

R_meanpH_model <- lmer(R ~ pH_mean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(R_meanpH_model)
summary(R_meanpH_model)

# R and range pH # 
R_rangepH_plot <- metadata %>%
  ggplot(aes(x = pH_rangemean, y = R)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in pH", y = "Respiration Rate") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 0.00, label.y = 0.07, size = 8) + 
  stat_cor(label.x = 0.00, label.y = 0.065, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
R_rangepH_plot
#ggsave(plot = R_rangepH_plot, filename = here("Output", "Biogeochem_Physio", "R_rangepH.png"), width = 14, height = 10)

R_rangepH_model <- lmer(R ~ pH_rangemean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(R_rangepH_model)
summary(R_rangepH_model)


# R and mean DOC # 
R_meanDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_mean, y = R)) + 
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
  stat_regline_equation(label.x = 150, label.y = 0.07, size = 8) + 
  stat_cor(label.x = 150, label.y = 0.065, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
R_meanDOC_plot
#ggsave(plot = R_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "R_meanDOC.png"), width = 14, height = 10)

R_meanDOC_model <- lmer(R ~ DOC_mean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(R_meanDOC_model)
summary(R_meanDOC_model)

# R and range DOC # 
R_rangeDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_rangemean, y = R)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in DOC", y = "Respiration Rate") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 0, label.y = 0.07, size = 8) + 
  stat_cor(label.x = 0, label.y = 0.065, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
R_rangeDOC_plot
#ggsave(plot = R_rangeDOC_plot, filename = here("Output", "Biogeochem_Physio", "R_rangeDOC.png"), width = 14, height = 10)

R_rangeDOC_model <- lmer(R ~ DOC_rangemean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(R_rangeDOC_model)
summary(R_rangeDOC_model)

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
GP_meanpH_plot <- metadata %>%
  ggplot(aes(x = pH_mean, y = GP)) + 
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
  stat_regline_equation(label.x = 8.045, label.y = 0.20, size = 5) + 
  stat_cor(label.x = 8.045, label.y = 0.19, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
GP_meanpH_plot
#ggsave(plot = GP_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "GP_meanpH.png"), width = 14, height = 12)

GP_meanpH_model <- lmer(GP ~ pH_mean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(GP_meanpH_model)
summary(GP_meanpH_model)

# GP and range pH # 
GP_rangepH_plot <- metadata %>%
  ggplot(aes(x = pH_rangemean, y = GP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in pH", y = "Gross Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 0.00, label.y = 0.20, size = 5) + 
  stat_cor(label.x = 0.00, label.y = 0.19, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
GP_rangepH_plot
#ggsave(plot = GP_rangepH_plot, filename = here("Output", "Biogeochem_Physio", "GP_rangepH.png"), width = 14, height = 12)

GP_rangepH_model <- lmer(GP ~ pH_rangemean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(GP_rangepH_model)
summary(GP_rangepH_model)

# GP and mean DOC # 
GP_meanDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_mean, y = GP)) + 
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
  stat_regline_equation(label.x = 150, label.y = 0.20, size = 5) + 
  stat_cor(label.x = 150, label.y = 0.19, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
GP_meanDOC_plot
#ggsave(plot = GP_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "GP_meanDOC.png"), width = 14, height = 12)

GP_meanDOC_model <- lmer(GP ~ DOC_mean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(GP_meanDOC_model)
summary(GP_meanDOC_model)

# GP and range DOC # 
GP_rangeDOC_plot <- metadata %>%
  ggplot(aes(x = DOC_rangemean, y = GP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in DOC", y = "Gross Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 40, label.y = 0.20, size = 5) + 
  stat_cor(label.x = 40, label.y = 0.19, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
GP_rangeDOC_plot
#ggsave(plot = GP_rangeDOC_plot, filename = here("Output", "Biogeochem_Physio", "GP_rangeDOC.png"), width = 14, height = 12)

GP_rangeDOC_model <- lmer(GP ~ DOC_rangemean + (1|CORAL_TANK_NUM) + (1|GENOTYPE), data = metadata)
check_model(GP_rangeDOC_model)
summary(GP_rangeDOC_model)

##### CALCULATE AVERAGE CORAL PHYSIO PARAMETERS PER TREATMENT #####

coral_physio_params_treatment <- metadata %>%
  group_by(TREATMENT) %>%
  summarize(avg_endos = mean(endo_per_cm2, na.rm = TRUE),
            endo_error = sd(endo_per_cm2, na.rm = TRUE)/sqrt(n()),
            avg_chl = mean(chla.ug.cm2, na.rm = TRUE),
            chl_error = sd(chla.ug.cm2, na.rm = TRUE)/sqrt(n()),
            avg_biomass = mean(mean_tissue_biomass_ng, na.rm = TRUE), 
            biomass_error = sd(mean_tissue_biomass_ng, na.rm = TRUE)/sqrt(n()),
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
