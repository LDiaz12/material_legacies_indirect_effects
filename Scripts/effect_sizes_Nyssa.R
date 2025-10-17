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
library(performance)

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
            mean_NEP = mean(NEP_mean, na.rm = TRUE),
            mean_chl_initial = mean(initial_chla, na.rm = TRUE))%>%
  left_join(raw_chem)

metadata$TREATMENT <- factor(metadata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

#metadata <- metadata %>%
#  left_join(raw_chem)

# standardize data to z scores first and log transform the chl, biomass, and endodata for normality
##- since they are all on different scales # 

metadata_scaled<-metadata %>%
  select(TREATMENT, GENOTYPE, chla_ug_cm2, endo_per_cm2, mean_tissue_biomass,GP, R, pH_mean, DOC_mean)%>%
  mutate_at(vars(chla_ug_cm2:mean_tissue_biomass), .funs = function(x){as.numeric(scale(log(x)))}) %>% # log transform then scale them
  mutate_at(vars(GP:DOC_mean), .funs = function(x){as.numeric(scale(x))}) # just scale these


# next we want to make a data sheet that pulls out the CI's and effect sizes # 
# make models for each parameter and save the estimate and CI summary table # 

## create lm of each mean physio parameter with grand mean pH ## 
## physiology ~ pH and DOC while controlling for genotype

endo_scale <- lmer(endo_per_cm2 ~ pH_mean +DOC_mean  +(1|GENOTYPE), data = metadata_scaled)
anova(endo_scale)
a <- coef(summary(endo_scale)) %>%
  as_tibble()%>%
  mutate(Chem = rownames(coef(summary(endo_scale))))%>%
  mutate(Parameter = "Endo_per_cm2")

a <-tibble(model_parameters(endo_scale)[1:3,])%>%
  mutate(Phys = "Endo_per_cm2")

chla_scale <- lmer(chla_ug_cm2 ~ pH_mean +DOC_mean  +(1|GENOTYPE), data = metadata_scaled)
b <- coef(summary(chla_scale)) %>%
  as_tibble()%>%
  mutate(Chem = rownames(coef(summary(chla_scale)))) %>%
mutate(Parameter = "Chla_ug_cm2")

b <- tibble(model_parameters(chla_scale)[1:3,])%>%
  mutate(Phys = "Chla_ug_cm2")

anova(chla_scale)

biomass_scale <- lmer(mean_tissue_biomass ~ pH_mean +DOC_mean  +(1|GENOTYPE), data = metadata_scaled)
c <- coef(summary(biomass_scale)) %>%
  as_tibble()%>%
  mutate(Chem = rownames(coef(summary(biomass_scale)))) %>%
  mutate(Parameter = "Mean_Tissue_Biomass")

c <- tibble(model_parameters(biomass_scale)[1:3,])%>%
  mutate(Phys = "Mean_Tissue_Biomass")
anova(biomass_scale)

#  breaks = c("Endo_per_cm2", "Chla_ug_cm2", "Mean_Tissue_Biomass", "R", "GP", "NP"),

R_scale <- lmer(R ~ pH_mean +DOC_mean  +(1|GENOTYPE), data = metadata_scaled)
d <- coef(summary(R_scale)) %>%
  as_tibble()%>%
  mutate(Chem = rownames(coef(summary(R_scale)))) %>%
  mutate(Parameter = "R")

d <- tibble(model_parameters(R_scale)[1:3,])%>%
  mutate(Phys = "R")

anova(R_scale)

GP_scale <- lmer(GP ~ pH_mean +DOC_mean  +(1|GENOTYPE), data = metadata_scaled)
e<-coef(summary(GP_scale)) %>%
  as_tibble()%>%
  mutate(Chem = rownames(coef(summary(GP_scale)))) %>%
  mutate(Parameter = "GP")

e <- tibble(model_parameters(GP_scale)[1:3,])%>%
  mutate(Phys = "GP")
anova(GP_scale)

## Now, rbind the summary data for each parameter into one data sheet ## 
chem.data1<-bind_rows(a,b,c,d,e) %>%
  filter(Parameter != "(Intercept)")%>% # remove intercept since we don't need it
  mutate(Parameter = ifelse(Parameter == "pH_mean","Mean pH","Mean DOC")) %>% # make prettier labels
  mutate(sig = ifelse(p< 0.055, 1, 0.5)) # create column with '1' for significant or '0' for not significant 


## Now we can make the plot of standardized effect sizes ## 
## effect size plot for grand mean pH and DOC separate, then combine using patchwork ##

chem.effect.plot1 <- ggplot(chem.data1, aes(x = Coefficient, y = Phys, col = Parameter, alpha=sig)) + 
  geom_point(size = 3, position = position_dodge(width = 0.4), shape = 19) + 
  labs(x = "Standardized Effect Size", y = "", 
       #title = "Mean pH", 
       color = ""
       ) +
  geom_errorbarh(aes(xmin = CI_low , xmax = CI_high), 
                 height = 0, position = position_dodge(0.4)) + 
  geom_vline(xintercept = 0, lty = 2) + 
  scale_alpha(range = c(0.4,1))+
  guides(alpha = "none")+
  scale_color_manual(values = c("#082a54","#800"))+
  scale_y_discrete(breaks = c("Endo_per_cm2", "Chla_ug_cm2", "Mean_Tissue_Biomass", "R", "GP"),
                   labels = c(expression(atop("Endosymbiont Density" , ~ (cells ~ 10^6 ~ cm^-2))), 
                              expression(atop("Chlorophyll-a Content" , ~ (µg ~ cm^-2))),
                              expression(atop("Mean Tissue Biomass" , ~ (mg ~ cm^-2))), 
                              expression(atop("Respiration Rate" , ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1))), 
                              expression(atop("Gross Photosynthesis" , ~ (µmol ~ O[2] ~ cm^-2 ~ hr^-1))))) +
 # scale_color_manual(values = c("Grand.Mean.pH" = "grey")) +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"), 
        axis.text.y = element_text(size = 11, color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 11)
        )
chem.effect.plot1


ggsave(plot = chem.effect.plot1, filename = here("Output", "effect_size_update.png"), width = 6, height = 5)

# here is where you get all your p values, coeffients, etc.
chem.data1
