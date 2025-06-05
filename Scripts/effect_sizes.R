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
  filter(!chla.ug.cm2 < -0.05,
         !mean_tissue_biomass > 0.00075) # points were greatly skewing data so was removed - also negative chl??

# standardize data to z scores first - since they are all on different scales # 
metadata$endo_per_cm2 <- scale(metadata$endo_per_cm2, center = TRUE, scale = TRUE)
metadata$chla.ug.cm2 <- scale(metadata$chla.ug.cm2, center = TRUE, scale = TRUE)
metadata$mean_tissue_biomass <- scale(metadata$mean_tissue_biomass, center = TRUE, scale = TRUE)
metadata$R <- scale(metadata$R, center = TRUE, scale = TRUE)
metadata$GP <- scale(metadata$GP, center = TRUE, scale = TRUE)
metadata$NP <- scale(metadata$NP, center = TRUE, scale = TRUE)

# next we want to make a data sheet that pulls out the CI's and effect sizes # 
# make models for each parameter and save the estimate and CI summary table # 
# FIRST WITH PH MEAN #
endo_scale <- lmer(endo_per_cm2 ~ pH_mean + (1|TANK_NUM), data = metadata)
a <- model_parameters(endo_scale)
anova(endo_scale)

chla_scale <- lmer(chla.ug.cm2 ~ pH_mean + (1|TANK_NUM), data = metadata)
b <- model_parameters(chla_scale)
anova(chla_scale)

biomass_scale <- lmer(mean_tissue_biomass ~ pH_mean + (1|TANK_NUM), data = metadata)
c <- model_parameters(biomass_scale)
anova(biomass_scale)

R_scale <- lmer(R ~ pH_mean + (1|TANK_NUM), data = metadata)
d <- model_parameters(R_scale)
anova(R_scale)

GP_scale <- lmer(GP ~ pH_mean + (1|TANK_NUM), data = metadata)
e <- model_parameters(GP_scale)
anova(GP_scale)

NP_scale <- lmer(NP ~ pH_mean + (1|TANK_NUM), data = metadata)
f <- model_parameters(NP_scale)
anova(NP_scale)

# THEN WITH DOC MEAN # 
endo_doc_scale <- lmer(endo_per_cm2 ~ DOC_mean + (1|TANK_NUM), data = metadata)
g <- model_parameters(endo_doc_scale)
anova(endo_doc_scale)

chla_doc_scale <- lmer(chla.ug.cm2 ~ DOC_mean + (1|TANK_NUM), data = metadata)
h <- model_parameters(chla_doc_scale)
anova(chla_doc_scale)

biomass_doc_scale <- lmer(mean_tissue_biomass ~ DOC_mean + (1|TANK_NUM), data = metadata)
i <- model_parameters(biomass_doc_scale)
anova(biomass_doc_scale)

R_doc_scale <- lmer(R ~ DOC_mean + (1|TANK_NUM), data = metadata)
j <- model_parameters(R_doc_scale)
anova(R_doc_scale)

GP_doc_scale <- lmer(GP ~ DOC_mean + (1|TANK_NUM), data = metadata)
k <- model_parameters(GP_doc_scale)
anova(GP_doc_scale)

NP_doc_scale <- lmer(NP ~ DOC_mean + (1|TANK_NUM), data = metadata)
l <- model_parameters(NP_doc_scale)
anova(NP_doc_scale)

## Now, rbind the summary data for each parameter into one data sheet ## 
chem.data <- rbind(a,b,c,d,e,f,g,h,i,j,k,l)
chem.data <- as.data.frame(chem.data)
# delete rows with intercept info - don't need 
chem.data <- chem.data[-c(1,3,4,5,7,8,9,11,12,13,15,16,17,19,20,21,23,24,25,27,28,29,31,32,33,35,36,37,39,40,41,43,44,45,47,48), ]

# create physio param column with labels
chem.data$phys.param <- c("Endo_per_cm2", "Chla_ug_cm2","Mean_Tissue_Biomass",
                        "R", "GP", "NP",
                        "Endo_per_cm2", "Chla_ug_cm2","Mean_Tissue_Biomass",
                        "R", "GP", "NP")
# create chem param column 
chem.data$chem.param <- c("Mean.pH", "Mean.pH", "Mean.pH", 
                        "Mean.pH", "Mean.pH", "Mean.pH",
                        "Mean.DOC", "Mean.DOC", "Mean.DOC", 
                        "Mean.DOC", "Mean.DOC", "Mean.DOC")
# create column with 'sig' or 'non' 
chem.data$sig <- c("non", "sig", "non", "non", "non", "non",
                   "non", "non", "non", "non", "non", "non")

# rename coefficient column to "estimate" 
colnames(chem.data)[2] <- "Estimate"

## Now we can make the plot of standardized effect sizes ## 
chem.effect.plot <- ggplot(chem.data, aes(x = Estimate, y = phys.param, col = chem.param, alpha = sig)) + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"), 
        axis.text.y = element_text(size = 11, color = "black")) + 
  geom_point(size = 3, alpha = 0.75) + 
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0) + 
  geom_vline(xintercept = 0, lty = 2) + 
  scale_y_discrete(breaks = c("Endo_per_cm2", "Chla_ug_cm2", "Mean_Tissue_Biomass", "R", "GP", "NP"),
                   labels = c("Endosymbiont Density \n(cells cm^-2)", "Chlorophyll a Content \n(ug cm^-2)",
                              "Mean Tissue Biomass \n(mg cm^-2)", "Respiration Rate", "Gross Photosynthesis",
                              "Net Photosynthesis")) +
  guides(alpha = FALSE) +
  scale_color_manual(labels = c("Mean.DOC", "Mean.pH"),
                       values = c("Mean.DOC" = "pink", "Mean.pH" = "lightblue")) +
  labs(col = "pH and DOC") + 
  xlab("Standardized Effect Sizes") + 
  ylab("")
chem.effect.plot
#ggsave(plot = chem.effect.plot, filename = here("Output", "effect_size_plot.png"), width = 9, height = 9)







