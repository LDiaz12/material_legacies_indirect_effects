library(tidyverse)
library(here)
library(agricolae)
library(lme4)
library(lmerTest)
library(sjPlot)
library(sjmisc)
library(effects)
library(sjstats)
library(effsize)
#install.packages("effsize")
library(parameters)
#install.packages("psycho")
library(psycho)
library(sjlabelled)
#install.packages("merDeriv")
library(merDeriv)

## read in metadata sheet # 
metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))
raw_chem <- read_csv(here("Data", "Chemistry", "Raw_Chem_Data.csv"))

metadata_raw_chem_means <- metadata %>%
  group_by(TANK_NUM, TREATMENT) %>%
  summarize(mean_endos = mean(endo_per_cm2, na.rm = TRUE),
            mean_chl = mean(chla_ug_cm2, na.rm = TRUE), 
            mean_biomass = mean(mean_tissue_biomass, na.rm = TRUE), 
            mean_GP = mean(GP, na.rm = TRUE),
            mean_R = mean(R, na.rm = TRUE),
            mean_pH = mean(pH_mean, na.rm = TRUE), 
            mean_DOC = mean(DOC_mean, na.rm = TRUE),
            deltapH = mean(deltapH_mean, na.rm = TRUE),
            deltaDOC = mean(deltaDOC_mean, na.rm = TRUE),
            mean_NEP = mean(NEP_mean, na.rm = TRUE)) %>%
  left_join(raw_chem)

metadata$TREATMENT <- factor(metadata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

metadata <- metadata %>%
  left_join(raw_chem)

# standardize data to z scores first - since they are all on different scales # 
metadata_raw_chem_means$mean_endos <- scale(metadata_raw_chem_means$mean_endos, center = TRUE, scale = TRUE)
metadata_raw_chem_means$mean_chl <- scale(metadata_raw_chem_means$mean_chl, center = TRUE, scale = TRUE)
metadata_raw_chem_means$mean_biomass <- scale(metadata_raw_chem_means$mean_biomass, center = TRUE, scale = TRUE)
metadata_raw_chem_means$mean_R <- scale(metadata_raw_chem_means$mean_R, center = TRUE, scale = TRUE)
metadata_raw_chem_means$mean_GP <- scale(metadata_raw_chem_means$mean_GP, center = TRUE, scale = TRUE)

# next we want to make a data sheet that pulls out the CI's and effect sizes # 
# make models for each parameter and save the estimate and CI summary table # 

## create lm of each mean physio parameter with grand mean pH ## 
endo_scale <- lm(mean_endos ~ grand_mean_pH, data = metadata_raw_chem_means)
a <- model_parameters(endo_scale)
anova(endo_scale)

chla_scale <- lm(mean_chl ~ grand_mean_pH, data = metadata_raw_chem_means)
b <- model_parameters(chla_scale)
anova(chla_scale)

biomass_scale <- lm(mean_biomass ~ grand_mean_pH, data = metadata_raw_chem_means)
c <- model_parameters(biomass_scale)
anova(biomass_scale)

R_scale <- lm(mean_R ~ grand_mean_pH, data = metadata_raw_chem_means)
d <- model_parameters(R_scale)
anova(R_scale)

GP_scale <- lm(mean_GP ~ grand_mean_pH, data = metadata_raw_chem_means)
e <- model_parameters(GP_scale)
anova(GP_scale)

## create lm of each mean physio parameter with grand mean DOC ## 
endo_doc_scale <- lm(mean_endos ~ grand_mean_DOC, data = metadata_raw_chem_means)
f <- model_parameters(endo_doc_scale)
anova(endo_doc_scale)

chla_doc_scale <- lm(mean_chl ~ grand_mean_DOC, data = metadata_raw_chem_means)
g <- model_parameters(chla_doc_scale)
anova(chla_doc_scale)

biomass_doc_scale <- lm(mean_biomass ~ grand_mean_DOC, data = metadata_raw_chem_means)
h <- model_parameters(biomass_doc_scale)
anova(biomass_doc_scale)

R_doc_scale <- lm(mean_R ~ grand_mean_DOC, data = metadata_raw_chem_means)
i <- model_parameters(R_doc_scale)
anova(R_doc_scale)

GP_doc_scale <- lm(mean_GP ~ grand_mean_DOC, data = metadata_raw_chem_means)
j <- model_parameters(GP_doc_scale)
anova(GP_doc_scale)

## Now, rbind the summary data for each parameter into one data sheet ## 
chem.data1 <- rbind(a,b,c,d,e) # grand mean pH 
chem.data2 <- rbind(f,g,h,i,j) # grand mean DOC

# make them a data frame 
chem.data1 <- as.data.frame(chem.data1)
chem.data2 <- as.data.frame(chem.data2)

# delete rows with intercept and obs info - don't need 
chem.data1 <- chem.data1[!grepl("Intercept|Observations", chem.data1$Parameter), (invert = TRUE), ]
chem.data2 <- chem.data2[!grepl("Intercept|Observations", chem.data2$Parameter), (invert = TRUE), ]

# create physio param column with labels
chem.data1$phys.param <- c("Endo_per_cm2", "Chla_ug_cm2","Mean_Tissue_Biomass",
                        "R", "GP")

chem.data2$phys.param <- c("Endo_per_cm2", "Chla_ug_cm2","Mean_Tissue_Biomass",
                           "R", "GP")

# create chem param column 
chem.data1$chem.param <- c("Grand.Mean.pH", "Grand.Mean.pH", "Grand.Mean.pH","Grand.Mean.pH", "Grand.Mean.pH")
                        
chem.data2$chem.param <- c("Grand.Mean.DOC", "Grand.Mean.DOC", "Grand.Mean.DOC", "Grand.Mean.DOC", "Grand.Mean.DOC")


# create column with '1' for significant or '0' for not significant 
chem.data1 <- chem.data1 %>%
  mutate(sig = ifelse(p < 0.05, "1", "0"))

chem.data2 <- chem.data2 %>%
  mutate(sig = ifelse(p < 0.05, "1", "0"))

# rename coefficient column to "estimate" 
colnames(chem.data1)[2] <- "Estimate"
colnames(chem.data2)[2] <- "Estimate"

## Now we can make the plot of standardized effect sizes ## 
## effect size plot for grand mean pH and DOC separate, then combine using patchwork ##

chem.effect.plot1 <- ggplot(chem.data1, aes(x = Estimate, y = phys.param, col = chem.param)) + 
  geom_point(size = 3) + 
  labs(x = "Standardized Effect Size", y = "", title = "Mean pH", color = "") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0) + 
  geom_vline(xintercept = 0, lty = 2) + 
  scale_y_discrete(breaks = c("Endo_per_cm2", "Chla_ug_cm2", "Mean_Tissue_Biomass", "R", "GP", "NP"),
                   labels = c(expression(atop("Endosymbiont Density" , ~ (cells ~ 10^6 ~ cm^-2))), 
                                         expression(atop("Chlorophyll-a Content" , ~ (µg ~ cm^-2))),
                              expression(atop("Mean Tissue Biomass" , ~ (mg ~ cm^-2))), 
                                         expression(atop("Respiration Rate" , ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1))), 
                              expression(atop("Gross Photosynthesis" , ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1))))) +
  scale_color_manual(values = c("Grand.Mean.pH" = "grey")) +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"), 
        axis.text.y = element_text(size = 11, color = "black"),
        legend.position = "none")
chem.effect.plot1
ggsave(plot = chem.effect.plot1, filename = here("Output", "effect_size_plot_pH.png"), width = 9, height = 9)

chem.effect.plot2 <- ggplot(chem.data2, aes(x = Estimate, y = phys.param, col = chem.param)) + 
  geom_point(size = 3) + 
  labs(x = "Standardized Effect Size", y = "", title = "Mean DOC", color = "") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0) + 
  geom_vline(xintercept = 0, lty = 2) + 
  scale_color_manual(values = c("Grand.Mean.DOC" = "grey")) +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"), 
        axis.text.y = element_blank(),
        legend.position = "none")
chem.effect.plot2
ggsave(plot = chem.effect.plot2, filename = here("Output", "effect_size_plot_DOC.png"), width = 9, height = 9)

effect_size_patch <- (chem.effect.plot1 + chem.effect.plot2) + plot_annotation(tag_levels = "a")
effect_size_patch
ggsave(plot = effect_size_patch, filename = here("Output", "effect_size_patch.png"), width = 9, height = 9)
