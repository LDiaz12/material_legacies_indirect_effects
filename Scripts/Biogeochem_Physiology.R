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

library(GGally)

## ANALYSIS OF BIOGEOCHEM IMPACTING CORAL PHYSIOLOGY ## 

# load in metadata sheet # 
metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))
metadata <- metadata %>%
  mutate(mean_tissue_biomass = ifelse(mean_tissue_biomass > 30, NA, mean_tissue_biomass))
raw_chem <- read_csv(here("Data", "Chemistry", "Raw_Chem_Data.csv"))

metadata_raw_chem_means <- metadata %>%
  group_by(TANK_NUM, TREATMENT) %>%
  summarize(mean_endos = mean(endo_per_cm2, na.rm = TRUE),
            se_endos = sd(endo_per_cm2, na.rm = TRUE)/sqrt(n()),
         mean_chl = mean(chla_ug_cm2, na.rm = TRUE), 
         se_chl = sd(chla_ug_cm2, na.rm = TRUE)/sqrt(n()),
         mean_biomass = mean(mean_tissue_biomass, na.rm = TRUE),
         se_biomass = sd(mean_tissue_biomass, na.rm = TRUE)/sqrt(n()),
         mean_GP = mean(GP, na.rm = TRUE),
         se_GP = sd(GP, na.rm = TRUE)/sqrt(n()),
         mean_R = mean(R, na.rm = TRUE),
         se_R = sd(R, na.rm = TRUE)/sqrt(n()),
         mean_pH = mean(pH_mean, na.rm = TRUE), 
         se_pH = sd(pH_mean, na.rm = TRUE)/sqrt(n()),
         mean_DOC = mean(DOC_mean, na.rm = TRUE),
         se_DOC = sd(DOC_mean, na.rm = TRUE)/sqrt(n()),
         deltapH = mean(deltapH_mean, na.rm = TRUE),
         deltaDOC = mean(deltaDOC_mean, na.rm = TRUE),
         mean_NEP = mean(NEP_mean, na.rm = TRUE)) %>%
  left_join(raw_chem) 
 # left_join(daily_min_max)

##last thing - check daily max and average 


####### ANOVAS OF PHYSIO PARAMS PER TREATMENT #######
metadata$TREATMENT <- factor(metadata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

# ENDOS #
endos_plot <- metadata %>% 
  ggplot(aes(x = TREATMENT, y = endo_per_cm2, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Endosymbiont Density" ~ (cells ~ x10^6 ~ cm^-2)))) +
  scale_x_discrete(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata, aes(x = TREATMENT, y = endo_per_cm2), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  coord_trans(y = "log")+
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
endos_plot
#ggsave(plot = endos_plot, filename = here("Output", "Biogeochem_Physio", "endos_plot.png"), width = 14, height = 10)

endos_model <- lmer(log(endo_per_cm2) ~ TREATMENT + (1|GENOTYPE), data= metadata)
check_model(endos_model)
summary(endos_model)
anova(endos_model) # non significant. p = 0.246

# CHL A #
chla_plot <- metadata %>% 
  ggplot(aes(x = TREATMENT, y = chla_ug_cm2, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Chlorophyll-a Content" ~ (µg ~ cm^-2)))) +
  scale_x_discrete(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata, aes(x = TREATMENT, y = chla_ug_cm2), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  coord_trans(y = "log")+
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
chla_plot
#ggsave(plot = chla_plot, filename = here("Output", "Biogeochem_Physio", "chla_plot.png"), width = 14, height = 10)

chla_model <- lmer(log(chla_ug_cm2) ~ TREATMENT + (1|GENOTYPE), data = metadata)
check_model(chla_model)
summary(chla_model)
anova(chla_model) # non significant. p = 0.259


# MEAN TISSUE BIOMASS # 
tissuebiomass_plot <- metadata %>% 
  ggplot(aes(x = TREATMENT, y = mean_tissue_biomass, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Tissue Biomass" ~ (mg ~ cm^-2)))) +
  scale_x_discrete(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = metadata, aes(x = TREATMENT, y = mean_tissue_biomass), alpha = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  coord_trans(y = "log")+
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
tissuebiomass_plot
#ggsave(plot = tissuebiomass_plot, filename = here("Output", "Biogeochem_Physio", "tissuebiomass_plot.png"), width = 14, height = 10)

tissuebiomass_model <- lmer(log(mean_tissue_biomass) ~ TREATMENT + (1|GENOTYPE), data = metadata)
check_model(tissuebiomass_model)
summary(tissuebiomass_model)
anova(tissuebiomass_model) #p = 0.147


# R #
R_plot <- metadata %>% 
  ggplot(aes(x = TREATMENT, y = R, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Respiration Rate" ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1)))) +
  scale_x_discrete(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched")) +
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

R_model <- lmer(R ~ TREATMENT + (1|GENOTYPE), data=metadata)
check_model(R_model)
summary(R_model)
anova(R_model) # p = 0.16


# GP # 
GP_plot <- metadata %>% 
  ggplot(aes(x = TREATMENT, y = GP, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Gross Photosynthesis" ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1)))) +
  scale_x_discrete(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched")) +
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

GP_model <- lmer(GP ~ TREATMENT + (1|GENOTYPE), data=metadata)
check_model(GP_model)
summary(GP_model)
anova(GP_model) # s. p = 0.61

physio_params_patch <- (endos_plot + chla_plot + tissuebiomass_plot)/(R_plot + GP_plot) + plot_annotation(tag_levels = "a")
physio_params_patch
ggsave(plot = physio_params_patch, filename = here("Output", "Supp_Fig_2.png"), width = 11, height = 11)

######## BIOGEOCHEM IMPACTS ON CHL-A ########
metadata_raw_chem_means$TREATMENT <- factor(metadata_raw_chem_means$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

# regression of chl-a and mean pH # 
chla_meanpH_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_pH, y = mean_chl)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x, color = "black") + 
  labs(x = expression("pH"[T]), y = expression("Chlorophyll a Content " ~ (µg ~ cm^-2))) +
  coord_trans(y = "log")+
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "bottom") +
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
                     values = c("blue", "darkgreen", "coral", "tan"))
chla_meanpH_plot
#ggsave(plot = chla_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "chla_meanpH.png"), width = 14, height = 10)

chla_meanpH_model <- lm(mean_chl ~ grand_mean_pH, data = metadata_raw_chem_means)
check_model(chla_meanpH_model)
anova(chla_meanpH_model)
summary(chla_meanpH_model)

# regression of chl-a and mean DOC # 
chla_meanDOC_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_DOC, y = mean_chl)) + 
  geom_point(aes(color = TREATMENT)) +  
  labs(x = expression("DOC ("~mu~"mol L"^-1~")"), y = expression("Chlorophyll a Content " ~ (µg ~ cm^-2))) +
  coord_trans(y = "log")+
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "bottom") + 
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
                     values = c("blue", "darkgreen", "coral", "tan"))
chla_meanDOC_plot
#ggsave(plot = chla_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "chla_meanDOC.png"), width = 14, height = 10)

chla_meanDOC_model <- lm(mean_chl ~ grand_mean_DOC, data = metadata_raw_chem_means)
check_model(chla_meanDOC_model)
anova(chla_meanDOC_model) 

######## BIOGEOCHEM IMPACTS ON ENDO DENSITY #######

# regression of endos and mean pH #
endo_meanpH_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_pH, y = mean_endos)) + 
  geom_point(aes(color = TREATMENT)) + 
  labs(x = expression("pH"[T]), y = expression("Endosymbiont Density " ~ (cells ~ 10^6 ~ cm^-2))) +
  coord_trans(y = "log")+
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"), 
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "bottom") + 
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
                     values = c("blue", "darkgreen", "coral", "tan"))
endo_meanpH_plot
#ggsave(plot = endo_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "endo_meanpH.png"), width = 14, height = 10)

endo_meanpH_model <- lm(mean_endos ~ grand_mean_pH, data = metadata_raw_chem_means)
check_model(endo_meanpH_model)
anova(endo_meanpH_model)

# regression of endos and mean DOC #
endo_meanDOC_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_DOC, y = mean_endos)) + 
  geom_point(aes(color = TREATMENT)) +  
  labs(x = expression("DOC ("~mu~"mol L"^-1~")"), y = expression("Endosymbiont Density " ~ (cells ~ 10^6 ~ cm^-2))) +
  coord_trans(y = "log")+
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "bottom") +
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
                     values = c("blue", "darkgreen", "coral", "tan"))
endo_meanDOC_plot
#ggsave(plot = endo_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "endo_meanDOC.png"), width = 14, height = 10)

endo_meanDOC_model <- lm(mean_endos ~ grand_mean_DOC, data = metadata_raw_chem_means)
check_model(endo_meanDOC_model)
anova(endo_meanDOC_model)

####### BIOGEOCHEM IMPACTS ON TISSUE BIOMASS ### #####

# regression of tissue biomass and mean pH # 
biomass_meanpH_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_pH, y = mean_biomass)) + 
  geom_point(aes(color = TREATMENT)) + 
  labs(x = expression("pH"[T]), y = expression("Mean Tissue Biomass " ~ (mg ~ cm^-2))) +
  coord_trans(y = "log")+
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "bottom") +  
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
                     values = c("blue", "darkgreen", "coral", "tan"))
biomass_meanpH_plot
#ggsave(plot = biomass_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "biomass_meanpH.png"), width = 14, height = 10)

biomass_meanpH_model <- lm(mean_biomass ~ grand_mean_pH, data = metadata_raw_chem_means)
check_model(biomass_meanpH_model)
anova(biomass_meanpH_model)

# regression of tissue biomass and mean DOC # 
biomass_meanDOC_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_DOC, y = mean_biomass)) + 
  geom_point(aes(color = TREATMENT)) + 
  coord_trans(y = "log")+
  labs(x = expression("DOC ("~mu~"mol L"^-1~")"), y = expression("Mean Tissue Biomass " ~ (mg ~ cm^-2))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "bottom") + 
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
                     values = c("blue", "darkgreen", "coral", "tan"))
biomass_meanDOC_plot
#ggsave(plot = biomass_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "biomass_meanDOC.png"), width = 14, height = 10)

biomass_meanDOC_model <- lm(mean_biomass ~ grand_mean_DOC, data = metadata_raw_chem_means)
check_model(biomass_meanDOC_model)
anova(biomass_meanDOC_model)

##### BIOGEOCHEM IMPACTS ON RESPIRATION RATE #####

# R and mean pH #
R_meanpH_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_pH, y = mean_R)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x, color = "black") + 
  labs(x = expression("pH"[T]), y = expression("Respiration Rate " ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "bottom") +  
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
                     values = c("blue", "darkgreen", "coral", "tan"))
R_meanpH_plot
#ggsave(plot = R_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "R_meanpH.png"), width = 14, height = 10)

R_meanpH_model <- lm(mean_R ~ grand_mean_pH, data = metadata_raw_chem_means)
check_model(R_meanpH_model)
summary(R_meanpH_model)
anova(R_meanpH_model)


# R and mean DOC # 
R_meanDOC_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_DOC, y = mean_R)) + 
  geom_point(aes(color = TREATMENT)) + 
  labs(x = expression("DOC ("~mu~"mol L"^-1~")"), y = expression("Respiration Rate " ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "bottom") + 
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
                     values = c("blue", "darkgreen", "coral", "tan"))
R_meanDOC_plot
#ggsave(plot = R_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "R_meanDOC.png"), width = 14, height = 10)

R_meanDOC_model <- lm(mean_R ~ grand_mean_DOC, data = metadata_raw_chem_means)
check_model(R_meanDOC_model)
anova(R_meanDOC_model)

##### BIOGEOCHEM IMPACTS ON GROSS PHOTOSYNTHESIS #####

# GP and mean pH # 
GP_meanpH_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_pH, y = mean_GP)) + 
  geom_point(aes(color = TREATMENT)) + 
  labs(x = expression("pH"[T]), y = expression("Gross Photosynthesis " ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "bottom") +  
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
                     values = c("blue", "darkgreen", "coral", "tan"))
GP_meanpH_plot
#ggsave(plot = GP_meanpH_plot, filename = here("Output", "Biogeochem_Physio", "GP_meanpH.png"), width = 14, height = 12)

GP_meanpH_model <- lm(mean_GP ~ grand_mean_pH, data = metadata_raw_chem_means)
check_model(GP_meanpH_model)
anova(GP_meanpH_model)

# GP and mean DOC # 
GP_meanDOC_plot <- metadata_raw_chem_means %>%
  ggplot(aes(x = grand_mean_DOC, y = mean_GP)) + 
  geom_point(aes(color = TREATMENT)) + 
  labs(x = expression("DOC ("~mu~"mol L"^-1~")"), y = expression("Gross Photosynthesis " ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),  
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = "bottom") +  
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched", "CCA-Enriched"),
                     values = c("blue", "darkgreen", "coral", "tan"))
GP_meanDOC_plot
#ggsave(plot = GP_meanDOC_plot, filename = here("Output", "Biogeochem_Physio", "GP_meanDOC.png"), width = 14, height = 12)

GP_meanDOC_model <- lm(mean_GP ~ grand_mean_DOC, data = metadata_raw_chem_means)
check_model(GP_meanDOC_model)
anova(GP_meanDOC_model)

meanpH_physio_patch <- (chla_meanpH_plot + endo_meanpH_plot + biomass_meanpH_plot)/(R_meanpH_plot + GP_meanpH_plot) + plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom") 
meanpH_physio_patch
ggsave(plot = meanpH_physio_patch, filename = here("Output", "Supp_Fig_3.png"), width = 15, height = 11)

meanDOC_physio_patch <- (chla_meanDOC_plot + endo_meanDOC_plot + biomass_meanDOC_plot)/(R_meanDOC_plot + GP_meanDOC_plot) + plot_annotation(tag_levels = "a") + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
meanDOC_physio_patch
ggsave(plot = meanDOC_physio_patch, filename = here("Output", "Supp_Fig_4.png"), width = 15, height = 11)


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

avg_GP_plot <- coral_physio_params_treatment %>%
  ggplot(aes(x = TREATMENT, y = avg_GP, color = TREATMENT)) + 
  geom_point() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_GP_plot

##### CALCULATE AVERAGE CORAL PHYSIO PARAMS PER TANK ##### 
coral_physio_params_tank <- metadata %>%
  group_by(TANK_NUM, TREATMENT) %>%
  summarize(avg_endos = mean(endo_per_cm2, na.rm = TRUE),
            avg_chl = mean(chla_ug_cm2, na.rm = TRUE), 
            avg_biomass = mean(mean_tissue_biomass, na.rm = TRUE), 
            avg_R = mean(R, na.rm = TRUE), 
            avg_GP = mean(GP, na.rm = TRUE),
            avg_NEC = mean(NEC_mean, na.rm = TRUE), 
            avg_NEP = mean(NEP_mean, na.rm = TRUE))
coral_physio_params_tank
write_csv(coral_physio_params_tank, here("Data", "Summary_Files", "coral_physio_params_tank.csv"))

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

avg_GP_plot <- coral_physio_params_tank %>%
  ggplot(aes(x = TREATMENT, y = avg_GP, color = TREATMENT)) + 
  geom_boxplot() +
  scale_color_manual(values = c("blue", "darkgreen", "coral", "tan"))
avg_GP_plot


## create correlation coefficient plot ### 

cor.test(metadata$GP, metadata$mean_tissue_biomass)

metadata_newnames <- metadata %>%
  rename(TissueBiomass = mean_tissue_biomass,
         Endosymbionts = endo_per_cm2,
         Chlorophyll = chla_ug_cm2)

physio <- metadata_newnames %>%
  select(TissueBiomass, R, GP, Endosymbionts, Chlorophyll)

# create custom function to only plot regression lines for significant p values # 
only_sig <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  cor_test <- cor.test(x, y)
  p_value <- cor_test$p.value
  
  p <- ggplot(data, mapping) + 
    geom_point(color = "coral", alpha = 0.5) 
  
  if(p_value < 0.05) {
    p <- p + geom_smooth(method = "lm", se = FALSE, color = "black", alpha = 0.3)
  }
  return(p)
}  

# need to create a custom function in order to get rid of the word "Corr" in the plot
custom_corr <- function(data, mapping, method, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y) 
  
  corr <- cor(x, y, method = method, use = "complete.obs")
  
  ggally_text(label = as.character(round(corr, 2)),
              mapping = aes(),
              xP = 0.5, vP = 0.5,
              color = "black",
              size = 8,
              ...
              )
}
corr_plot <- ggpairs(data = physio,
                     upper = list(continuous = wrap(custom_corr, method = "pearson")),
                     lower = list(continuous = wrap(only_sig)),
                     diag = list(continuous = wrap("densityDiag", fill = "lightblue"))) + 
  theme_minimal()
corr_plot <- corr_plot + 
  theme(axis.text.x = element_text(size = 12),
         axis.text.y = element_text(size = 12))
corr_plot
ggsave(plot = corr_plot, filename = here("Output", "physio_corr_coeff_plot.png"))

#this also creates the plot but with the "Corr" kept in 
corr_coeff_plot <- ggpairs(data = physio,
        upper = list(continuous = wrap("cor", method = "pearson", size = 6, color = "black")), 
        lower = list(continuous = wrap("smooth", color = "coral"), alpha = 0.3),
        diag = list(continuous = wrap("densityDiag", fill = "lightblue"))) +
  theme_minimal()

corr_coeff_plot

#ggsave(plot = corr_coeff_plot, filename = here("Output", "physio_corr_coeff_plot.png"))
