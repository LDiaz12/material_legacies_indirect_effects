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

metadata$TREATMENT <- factor(metadata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

# standardize data to z scores first - since they are all on different scales # 
metadata$endo_per_cm2 <- scale(metadata$endo_per_cm2, center = TRUE, scale = TRUE)
metadata$chla.ug.cm2 <- scale(metadata$chla.ug.cm2, center = TRUE, scale = TRUE)
metadata$mean_tissue_biomass_ng <- scale(metadata$mean_tissue_biomass_ng, center = TRUE, scale = TRUE)
metadata$R <- scale(metadata$R, center = TRUE, scale = TRUE)
metadata$GP <- scale(metadata$GP, center = TRUE, scale = TRUE)

# next we want to make a data sheet that pulls out the CI's and effect sizes # 
# make models for each parameter and save the estimate and CI summary table # 

# PH MEAN #
endo_scale <- lmer(endo_per_cm2 ~ pH_mean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
a <- model_parameters(endo_scale)
anova(endo_scale)

chla_scale <- lmer(chla.ug.cm2 ~ pH_mean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
b <- model_parameters(chla_scale)
anova(chla_scale)

biomass_scale <- lmer(mean_tissue_biomass_ng ~ pH_mean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
c <- model_parameters(biomass_scale)
anova(biomass_scale)

R_scale <- lmer(R ~ pH_mean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
d <- model_parameters(R_scale)
anova(R_scale)

GP_scale <- lmer(GP ~ pH_mean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
e <- model_parameters(GP_scale)
anova(GP_scale)

# DOC MEAN # 
endo_doc_scale <- lmer(endo_per_cm2 ~ DOC_mean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
f <- model_parameters(endo_doc_scale)
anova(endo_doc_scale)

chla_doc_scale <- lmer(chla.ug.cm2 ~ DOC_mean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
g <- model_parameters(chla_doc_scale)
anova(chla_doc_scale)

biomass_doc_scale <- lmer(mean_tissue_biomass_ng ~ DOC_mean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
h <- model_parameters(biomass_doc_scale)
anova(biomass_doc_scale)

R_doc_scale <- lmer(R ~ DOC_mean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
i <- model_parameters(R_doc_scale)
anova(R_doc_scale)

GP_doc_scale <- lmer(GP ~ DOC_mean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
j <- model_parameters(GP_doc_scale)
anova(GP_doc_scale)

# PH RANGE # 
endo_scale2 <- lmer(endo_per_cm2 ~ pH_rangemean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
k <- model_parameters(endo_scale2)
anova(endo_scale2)

chla_scale2 <- lmer(chla.ug.cm2 ~ pH_rangemean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
l <- model_parameters(chla_scale2)
anova(chla_scale2)

biomass_scale2 <- lmer(mean_tissue_biomass_ng ~ pH_rangemean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
m <- model_parameters(biomass_scale2)
anova(biomass_scale2)

R_scale2 <- lmer(R ~ pH_rangemean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
n <- model_parameters(R_scale2)
anova(R_scale2)

GP_scale2 <- lmer(GP ~ pH_rangemean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
o <- model_parameters(GP_scale2)
anova(GP_scale2)

# DOC RANGE # 
endo_doc_scale2 <- lmer(endo_per_cm2 ~ DOC_rangemean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
p <- model_parameters(endo_doc_scale2)
anova(endo_doc_scale2)

chla_doc_scale2 <- lmer(chla.ug.cm2 ~ DOC_rangemean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
q <- model_parameters(chla_doc_scale2)
anova(chla_doc_scale2)

biomass_doc_scale2 <- lmer(mean_tissue_biomass_ng ~ DOC_rangemean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
r <- model_parameters(biomass_doc_scale2)
anova(biomass_doc_scale2)

R_doc_scale2 <- lmer(R ~ DOC_rangemean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
s <- model_parameters(R_doc_scale2)
anova(R_doc_scale2)

GP_doc_scale2 <- lmer(GP ~ DOC_rangemean + (1|TANK_NUM) + (1|GENOTYPE), data = metadata)
t <- model_parameters(GP_doc_scale2)
anova(GP_doc_scale2)

## Now, rbind the summary data for each parameter into one data sheet ## 
chem.data1 <- rbind(a,b,c,d,e,k,l,m,n,o)
chem.data2 <- rbind(f,g,h,i,j,p,q,r,s,t)

chem.data1 <- as.data.frame(chem.data1)
chem.data2 <- as.data.frame(chem.data2)

# delete rows with intercept and obs info - don't need 
chem.data1 <- chem.data1[!grepl("Intercept|Observations", chem.data1$Parameter), (invert = TRUE), ]
chem.data2 <- chem.data2[!grepl("Intercept|Observations", chem.data2$Parameter), (invert = TRUE), ]

# create physio param column with labels
chem.data1$phys.param <- c("Endo_per_cm2", "Chla_ug_cm2","Mean_Tissue_Biomass",
                        "R", "GP",
                        "Endo_per_cm2", "Chla_ug_cm2","Mean_Tissue_Biomass",
                        "R", "GP")

chem.data2$phys.param <- c("Endo_per_cm2", "Chla_ug_cm2","Mean_Tissue_Biomass",
                           "R", "GP",
                           "Endo_per_cm2", "Chla_ug_cm2","Mean_Tissue_Biomass",
                           "R", "GP")

# create chem param column 
chem.data1$chem.param <- c("pH.Mean", "pH.Mean", "pH.Mean","pH.Mean", "pH.Mean",
                           "pH.Range", "pH.Range", "pH.Range", "pH.Range", "pH.Range")
                        
chem.data2$chem.param <- c("DOC.Mean", "DOC.Mean", "DOC.Mean", "DOC.Mean", "DOC.Mean",
                           "DOC.Range", "DOC.Range", "DOC.Range", "DOC.Range", "DOC.Range")


# create column with '1' for significant or '0' for not significant 
chem.data1 <- chem.data1 %>%
  mutate(sig = ifelse(p < 0.05, "1", "0"))

chem.data2 <- chem.data2 %>%
  mutate(sig = ifelse(p < 0.05, "1", "0"))

# rename coefficient column to "estimate" 
colnames(chem.data1)[2] <- "Estimate"
colnames(chem.data2)[2] <- "Estimate"

## Now we can make the plot of standardized effect sizes ## 
## effect size plot for mean and range separate then combine with patchwork ##

chem.effect.plot1 <- ggplot(chem.data1, aes(x = Estimate, y = phys.param, col = chem.param)) + 
  geom_point(size = 3) + 
  labs(x = "Standardized Effect Size", y = "", title = "pH", color = "") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0) + 
  geom_vline(xintercept = 0, lty = 2) + 
  scale_y_discrete(breaks = c("Endo_per_cm2", "Chla_ug_cm2", "Mean_Tissue_Biomass", "R", "GP", "NP"),
                   labels = c(expression(atop("Endosymbiont Density" , ~ (cells ~ cm^-2))), 
                                         expression(atop("Chlorophyll-a Content" , ~ (µg ~ cm^-2))),
                              expression(atop("Mean Tissue Biomass" , ~ (ng ~ cm^-2))), 
                                         expression(atop("Respiration Rate" , ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1))), 
                              expression(atop("Gross Photosynthesis" , ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1))))) +
  scale_color_manual(values = c("pH.Mean" = "black", "pH.Range" = "grey")) +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"), 
        axis.text.y = element_text(size = 11, color = "black"),
        legend.position = "none")
chem.effect.plot1
#ggsave(plot = chem.effect.plot1, filename = here("Output", "effect_size_plot_pH.png"), width = 9, height = 9)

chem.effect.plot2 <- ggplot(chem.data2, aes(x = Estimate, y = phys.param, col = chem.param)) + 
  geom_point(size = 3) + 
  labs(x = "Standardized Effect Size", y = "", title = "DOC", color = "") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0) + 
  geom_vline(xintercept = 0, lty = 2) + 
  scale_color_manual(values = c("DOC.Mean" = "black", "DOC.Range" = "grey")) +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"), 
        axis.text.y = element_blank(),
        legend.position = "none")
chem.effect.plot2
#ggsave(plot = chem.effect.plot2, filename = here("Output", "effect_size_plot_DOC.png"), width = 9, height = 9)

effect_size_patch <- (chem.effect.plot1 + chem.effect.plot2) + plot_annotation(tag_levels = "a")
effect_size_patch
#ggsave(plot = effect_size_patch, filename = here("Output", "effect_size_patch.png"), width = 9, height = 9)
