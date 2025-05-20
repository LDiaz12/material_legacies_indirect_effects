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

## ANALYSIS OF BIOGEOCHEM IMPACTING CORAL PHYSIOLOGY ## 

# load in metadata sheet # 
metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))

### BIOGEOCHEM IMPACTS ON CHL-A ###
metadata <- metadata %>%
  filter(!chla.ug.cm2 < -0.05,
         !mean_tissue_biomass > 0.00075) # points were greatly skewing data so was removed - also negative chl??

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

chla_meanpH_model <- lm(chla.ug.cm2 ~ pH_mean, data = metadata)
check_model(chla_meanpH_model)
summary(chla_meanpH_model) # barely significant at p = 0.0202

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

chla_rangepH_model <- lm(chla.ug.cm2 ~ pH_rangemean, data = metadata)
check_model(chla_rangepH_model)
summary(chla_rangepH_model) # significant at p = 0.04

# regression of chl-a and mean TA # 
chla_meanTA_plot <- metadata %>%
  ggplot(aes(x = TA_mean, y = chla.ug.cm2)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  stat_regline_equation(label.x = 2320, label.y = 0.14, size = 8) + 
  stat_cor(label.x = 2320, label.y = 0.13, size = 8) +
  labs(x = "Daily Mean TA", y = "Chlorophyll-a ug/cm2") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
chla_meanTA_plot
#ggsave(plot = chla_meanTA_plot, filename = here("Output", "Biogeochem_Physio", "chla_meanTA.png"), width = 14, height = 10)

# regression of chl-a and mean range in TA # 
chla_rangeTA_plot <- metadata %>%
  ggplot(aes(x = TA_rangemean, y = chla.ug.cm2)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in TA", y = "Chlorophyll-a ug/cm2") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = -35, label.y = 0.14, size = 8) + 
  stat_cor(label.x = -35, label.y = 0.13, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
chla_rangeTA_plot
#ggsave(plot = chla_rangeTA_plot, filename = here("Output", "Biogeochem_Physio", "chla_rangeTA.png"), width = 14, height = 10)

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


### BIOGEOCHEM IMPACTS ON ENDO DENSITY ###

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

# regression of endos and mean TA # 
endo_meanTA_plot <- metadata %>%
  ggplot(aes(x = TA_mean, y = endo_per_cm2)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean TA", y = "Endosymbiont Density (counts/cm2)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 2280, label.y = 1.3, size = 8) + 
  stat_cor(label.x = 2280, label.y = 1.2, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
endo_meanTA_plot
#ggsave(plot = endo_meanTA_plot, filename = here("Output", "Biogeochem_Physio", "endo_meanTA.png"), width = 14, height = 10)

# regression of endos and daily mean range in TA # 
endo_rangeTA_plot <- metadata %>%
  ggplot(aes(x = TA_rangemean, y = endo_per_cm2)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in TA", y = "Endosymbiont Density (counts/cm2)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = -75, label.y = 1.3, size = 8) + 
  stat_cor(label.x = -75, label.y = 1.2, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
endo_rangeTA_plot
#ggsave(plot = endo_rangeTA_plot, filename = here("Output", "Biogeochem_Physio", "endo_rangeTA.png"), width = 14, height = 10)

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

endo_rangeDOC_model <- lm(endo_per_cm2 ~ DOC_rangemean, data = metadata)
check_model(endo_rangeDOC_model)
summary(endo_rangeDOC_model) # significant at p = 0.023


### BIOGEOCHEM IMPACTS ON TISSUE BIOMASS ### 

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

# regression of tissue biomass and mean TA # 
biomass_meanTA_plot <- metadata %>%
  ggplot(aes(x = TA_mean, y = mean_tissue_biomass)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean TA", y = "Mean Tissue Biomass") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 2280, label.y = 4e-04, size = 8) + 
  stat_cor(label.x = 2280, label.y = 3.5e-04, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
biomass_meanTA_plot
#ggsave(plot = biomass_meanTA_plot, filename = here("Output", "Biogeochem_Physio", "biomass_meanTA.png"), width = 14, height = 10)

# regression of tissue biomass and range in TA # 
biomass_rangeTA_plot <- metadata %>%
  ggplot(aes(x = TA_rangemean, y = mean_tissue_biomass)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in TA", y = "Mean Tissue Biomass") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = -80, label.y = 4e-04, size = 8) + 
  stat_cor(label.x = -80, label.y = 3.5e-04, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
biomass_rangeTA_plot
#ggsave(plot = biomass_rangeTA_plot, filename = here("Output", "Biogeochem_Physio", "biomass_rangeTA.png"), width = 14, height = 10)

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

# R and mean TA # 
R_meanTA_plot <- metadata %>%
  ggplot(aes(x = TA_mean, y = R)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean TA", y = "Respiration Rate") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 2280, label.y = 0.07, size = 8) + 
  stat_cor(label.x = 2280, label.y = 0.065, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
R_meanTA_plot
#ggsave(plot = R_meanTA_plot, filename = here("Output", "Biogeochem_Physio", "R_meanTA.png"), width = 14, height = 10)

# R and range TA # 
R_rangeTA_plot <- metadata %>%
  ggplot(aes(x = TA_rangemean, y = R)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in TA", y = "Respiration Rate") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = -75, label.y = 0.07, size = 8) + 
  stat_cor(label.x = -75, label.y = 0.065, size = 8) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
R_rangeTA_plot
#ggsave(plot = R_rangeTA_plot, filename = here("Output", "Biogeochem_Physio", "R_rangeTA.png"), width = 14, height = 10)

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

# NP and mean TA # 
NP_meanTA_plot <- metadata %>%
  ggplot(aes(x = TA_mean, y = NP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean TA", y = "Net Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 2280, label.y = 0.14, size = 5) + 
  stat_cor(label.x = 2280, label.y = 0.13, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
NP_meanTA_plot
#ggsave(plot = NP_meanTA_plot, filename = here("Output", "Biogeochem_Physio", "NP_meanTA.png"), width = 14, height = 12)

# NP and range TA # 
NP_rangeTA_plot <- metadata %>%
  ggplot(aes(x = TA_rangemean, y = NP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in TA", y = "Net Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = -75, label.y = 0.14, size = 5) + 
  stat_cor(label.x = -75, label.y = 0.13, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
NP_rangeTA_plot
#ggsave(plot = NP_rangeTA_plot, filename = here("Output", "Biogeochem_Physio", "NP_rangeTA.png"), width = 14, height = 12)

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

# GP and mean TA #
GP_meanTA_plot <- metadata %>%
  ggplot(aes(x = TA_mean, y = GP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean TA", y = "Gross Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 2280, label.y = 0.20, size = 5) + 
  stat_cor(label.x = 2280, label.y = 0.19, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
GP_meanTA_plot
#ggsave(plot = GP_meanTA_plot, filename = here("Output", "Biogeochem_Physio", "GP_meanTA.png"), width = 14, height = 12)

# GP and range TA # 
GP_rangeTA_plot <- metadata %>%
  ggplot(aes(x = TA_rangemean, y = GP)) + 
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) + 
  labs(x = "Daily Mean Range in TA", y = "Gross Photosynthesis") +
  theme(axis.text.x = element_text(size = 15, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = -75, label.y = 0.20, size = 5) + 
  stat_cor(label.x = -75, label.y = 0.19, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
GP_rangeTA_plot
#ggsave(plot = GP_rangeTA_plot, filename = here("Output", "Biogeochem_Physio", "GP_rangeTA.png"), width = 14, height = 12)

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

##### CALCULATE AVERAGE CORAL PHYSIO PARAMETERS PER TREATMENT #####

coral_physio_params_treatment <- metadata %>%
  group_by(TREATMENT) %>%
  summarize(avg_endos = mean(endo_per_cm2, na.rm = TRUE),
            avg_chl = mean(chla.ug.cm2, na.rm = TRUE), 
            avg_biomass = mean(mean_tissue_biomass, na.rm = TRUE), 
            avg_R = mean(R, na.rm = TRUE), 
            avg_NP = mean(NP, na.rm = TRUE), 
            avg_GP = mean(GP, na.rm = TRUE))

coral_physio_params_treatment$TREATMENT <- factor(coral_physio_params_treatment$TREATMENT, levels = c("Control", "Algae_Dom",
                                                                                                      "Coral_Dom", "Rubble_Dom"))
avg_endos_plot <- coral_physio_params_treatment %>%
  ggplot(aes(x = TREATMENT, y = avg_endos, color = TREATMENT)) + 
  geom_point() + 
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
            avg_GP = mean(GP, na.rm = TRUE))

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
