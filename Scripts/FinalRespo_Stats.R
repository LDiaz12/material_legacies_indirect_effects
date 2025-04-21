library(tidyverse)
library(here)
library(performance)
library(emmeans)
library(ggrepel)

final_respo_norm <- read_csv(here("Data", "RespoFiles", "Final", "RespoR_Normalized_FinalRates.csv"))
metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))
clean_pH_data <- read_csv(here("Data", "Chemistry", "Cleaned_pH_Data_per_Treatment.csv"))
clean_pH_data_FULL <- read_csv(here("Data", "Chemistry", "Cleaned_pH_Data_FULL.csv"))
full_chem_data <- read_csv(here("Data", "Chemistry", "Full_Chem_Data.csv"))
pH_day_night_means <- read_csv(here("Data", "Chemistry", "pH_day_night_means.csv"))
physio_metadata <- read_csv(here("Data", "MO24BEAST_physio_metadata.csv"))

####################################
## RESPO AND PH DATA JOIN ##
pH_summary <- clean_pH_data_FULL %>%
  group_by(TREATMENT, TANK_NUM) %>%
  summarise(pH_mean = mean(pH_dailymean, na.rm = TRUE), 
            pH_range = mean(pH_range, na.rm = TRUE))


metadata$TREATMENT <- factor(metadata$TREATMENT, levels = c("Control", "Coral_Dom", "Algae_Dom", "Rubble_Dom"))

respo_meta_pHmeans <- final_respo_norm %>%
  right_join(metadata) %>% 
  select(TANK_NUM, GENOTYPE, TREATMENT, GP, R, NP) %>%
  left_join(pH_day_night_means)

#####################################
## GROSS PHOTOSYNTHETIC RATE ##
final_respo_GP <- final_respo_meta_join %>%
  ggplot(aes(x=TREATMENT, y = GP, color = TREATMENT)) +
  labs(x = "Community Tank", y = "Coral Gross Photosynthetic Rate (umol/cm2/hr)") +
  geom_jitter(data = final_respo_meta_join, aes(x = TREATMENT, y = GP), alpha = 0.7) +
  geom_text_repel(aes(label = CORAL_NUM)) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_respo_GP
# from here, I can see that the Control and Coral Dominated tanks have the lowest coral GP rates from the corals in the response tanks
# Algae dominated communities have the widest spread of respiration rates, from 0.05 - 0.22; could be due to the highest mean pH and 
# range in pH seen throughout the experiment -- BUT WHY? Stress? More acidic conditions = more stress = lower GP? Does higher pH = higher GP?

final_pHmean_GP <- final_respo_meta_join %>%
  ggplot(aes(x=pH_mean, y = GP, color = TREATMENT)) + 
  labs(x = "pH Mean", y = "Gross Photosynthesis (umol/cm2/hr)") +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHmean_GP # lots of warnings from the geom_smooth... 

# general trends: looks like the algae dominated treatments had a higher mean pH throughout the experiment (save for 
# a few outliers); with their corals experiencing higher GP when compared to Controls. 
# Control treatments stayed fairly consistent around 8.04 - 8.045 (with a couple outliers); Control corals fairly consistent 
# GP between 0.05 - 0.10 (with a few outliers). 
# Coral dominated communities seem to have a spread of mean pH throughout the tanks (tank effect?); their corals have GP rates 
# typically between 0.05 - 0.15. 
# Rubble dominated communities have the greatest range of pH values (~8.03. 8.04, 8.055, and 8.06); these corals however
# seem to stay between 0.075 - 0.15 

# final GP and mean pH from DAY sampling time # 
final_pHday_GP <- respo_meta_pHmeans %>%
  ggplot(aes(x=pH_day_mean, y = GP, color = TREATMENT)) + 
  labs(x = "Day pH Mean", y = "Gross Photosynthesis (umol/cm2/hr)") +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHday_GP

# final GP and mean pH from NIGHT sampling time # 
final_pHnight_GP <- respo_meta_pHmeans %>%
  ggplot(aes(x=pH_night_mean, y = GP, color = TREATMENT)) + 
  labs(x = "Night pH Mean", y = "Gross Photosynthesis (umol/cm2/hr)") +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHnight_GP

# final respo and pH daily range and GP 
final_pHrange_GP <- final_respo_meta_join %>%
  ggplot(aes(x=pH_range, y = GP, color = TREATMENT)) + # color = treatment
  labs(x = "pH Range", y = "Gross Photosynthesis (umol/cm2/hr)") +
  geom_point() + 
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHrange_GP

# here I'm seeing more of a distinct pattern compared to the mean pH: Control tanks had the least difference in pH range
# throughout sampling days (~0); corals in these treatments had GP rates between 0.05 - 0.15. 
# Coral dom communities had pH range between 0.065 - 0.125, and corals had GP rates between 0.05-0.15 like controls.
# Rubble dom had pH range between 0.08 - 0.16 and coral GP rates between 0.075 - 0.15. 
# Algae dom communities experienced the highest range of pH throughout the experimental period, reaching a pH 
# difference from 0.10 up to 0.175; coral GP rates were between 0.05, all the way up to almost 0.175. 

##################################
## RESPIRATION RATE ##
final_respo_R <- final_respo_meta_join %>%
  ggplot(aes(x=TREATMENT, y = R, color = TREATMENT)) +
  labs(x = "Community Tank", y = "Coral Respiration Rate (umol/cm2/hr)") + 
  geom_jitter(data = final_respo_meta_join, aes(x = TREATMENT, y = GP), alpha = 0.7) +
  geom_text_repel(aes(label = CORAL_NUM)) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_respo_R
# similar trends as seen with GP, again with Algae dominated communities having a larger spread of respiration rates

final_pHmean_R <- final_respo_meta_join %>%
  ggplot(aes(x=pH_mean, y = R, color = TREATMENT)) + 
  labs(x = "pH Mean", y = "Respiration Rate (umol/cm2/hr)") +
  geom_point() + 
  geom_smooth(method = "lm", formula = y~x) +
  #facet_wrap(~GENOTYPE) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHmean_R # more geom_smooth warnings... 

# looks like algae and rubble dominated communities had higher mean pH and higher respiration rates (with a lot more 
# spread of respiration rates)

# final R and DAY pH # 

final_pHday_R <- respo_meta_pHmeans %>%
  ggplot(aes(x=pH_day_mean, y = R, color = TREATMENT)) + 
  labs(x = "Day pH Mean", y = "Gross Photosynthesis (umol/cm2/hr)") +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHday_R

# final R and NIGHT pH # 
final_pHnight_R <- respo_meta_pHmeans %>%
  ggplot(aes(x=pH_night_mean, y = R, color = TREATMENT)) + 
  labs(x = "Night pH Mean", y = "Gross Photosynthesis (umol/cm2/hr)") +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHnight_R


# final R and pH daily range #
final_pHrange_R <- final_respo_meta_join %>%
  ggplot(aes(x=pH_range, y = R, color = TREATMENT)) + 
  labs(x = "pH Range", y = "Respiration Rate (umol/cm2/hr)") +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHrange_R
# more of a clear delineation here of the control community having the lowest pH range, similar to GP plots (~0)
# again, algae dom communities having a wide spread of respiration rates from 0.02 - 0.07
# can see some pretty clear groupings among treatments
# lower range in pH typically = lower respiration rate, whereas higher range in pH = higher respiration rate 


####################################
## NET PHOTOSYNTHESIS ## 
final_respo_NP <- final_respo_meta_join %>%
  ggplot(aes(x=TREATMENT, y = NP, color = TREATMENT)) +
  labs(x = "Community Tank", y = "Net Photosynthesis (umol/cm2/hr)") +
  geom_jitter(data = final_respo_meta_join, aes(x = TREATMENT, y = NP), alpha = 0.7) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_respo_NP
# similar trends to GP and R: control tanks had the lowest NP, followed by coral dom, rubble dom, then algae dom 
# with the highest rates of NP
# coral dom communities seem to have the least amount of variance in rates, while algae dom again has a wider spread

final_pHmean_NP <- final_respo_meta_join %>%
  ggplot(aes(x=pH_mean, y = NP, color = TREATMENT)) + 
  labs(x = "pH Mean", y = "Net Photosynthesis (umol/cm2/hr)") +
  geom_point() + 
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHmean_NP
# algae and rubble dom among the communities with the highest pH values, and corals from these tanks have some of 
# the highest NP rates
# coral dom communities seem to have NP values between 0.045 - 0.10

final_pHrange_NP <- final_respo_meta_join %>%
  ggplot(aes(x=pH_range, y = NP, color = TREATMENT)) + 
  labs(x = "pH Range", y = "Net Photosynthesis (umol/cm2/hr)") +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHrange_NP
# similar trends as seen before with control tanks; algae dom communities have the highest variance in NP while coral 
# dominated communities seem to have the least variance 

####################################
## GP/R: Ratio to assess coral energy production/usage. High GP/R = producing more energy; low = using more energy ##
final_respo_GP_R <- final_respo_meta_join %>%
  ggplot(aes(x=TREATMENT, y = GP/R, color = TREATMENT)) + 
  labs(x = "Community Tank", y = "Gross Photosynthesis/Respiration Ratio") +
  geom_jitter(data = final_respo_meta_join, aes(x = TREATMENT, y = GP/R), alpha = 0.7) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_respo_GP_R
# such a spread of data here: two of the rubble dom communities look like they may be clear outliers, but not sure
# coral dom communities look pretty consistent (not super high variance) around ratios of 2 -4 (is this typical for P rus?)
# maybe one outlier for algae dom community which could be skewing the data/plot 
# lots of variance with the control tanks which is really interesting... 

final_pHmean_GP_R <- final_respo_meta_join %>%
  ggplot(aes(x=pH_mean, y = GP/R, color = TREATMENT)) + 
  labs(x = "pH Mean", y = "Gross Photosynthesis/Respiration Ratio") +
  geom_point() + 
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHmean_GP_R
# data kinda all over the place, but can also see some trends here
# algae and rubble dom communities have higher mean pH and GP/R of 2 - 4
# MOST coral and some rubble dom communities have lower mean pH, and also have GP/R of 2 - 4 
# a few outliers here which were seen in the previous plot, so may be skewing 

final_pHrange_GP_R <- final_respo_meta_join %>%
  ggplot(aes(x=pH_range, y = GP/R, color = TREATMENT)) + 
  labs(x = "pH Range", y = "Gross Photosynthesis/Respiration Ratio") +
  geom_point() + 
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
final_pHrange_GP_R
# more clear delineation of groupings as seen with the other pH range plots 
# control tanks with ~0 range in pH values, but lots of variance for GP/R 
# coral dom communities around 0.06 - 0.12 range in pH values and GP/R ratios of 2 - 4 
# both rubble and algae dom communities between 0.10 - 0.17 range in pH values, and similarly have GP/R ratios of 2 - 4

R_GP_plot <- final_respo_meta_join %>%
  ggplot(aes(x = R, y = GP, color = TREATMENT)) + 
  geom_point() + 
  geom_smooth(method = "lm", formula = y~x)
R_GP_plot


#plot NEP per treatment# 
physio_metadata$TREATMENT <- factor(physio_metadata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

treatment_NEP_plot <- physio_metadata %>%
  ggplot(aes(x=TREATMENT, y = NEP, color = TREATMENT)) +
  labs(x = "Community Tank", y = "Net Ecosystem Production (mmol C)") +
  geom_jitter(data = physio_metadata, aes(x = TREATMENT, y = NEP), alpha = 0.7) +
  #geom_text_repel(aes(label = CORAL_NUM)) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
treatment_NEP_plot
#ggsave(plot = treatment_NEP_plot, filename = here("Output", "NEP_treatment_plot.png"), width = 9, height = 6)

# calculate mean NEP per treatment #
NEP_sum_data <- physio_metadata %>% 
  group_by(TREATMENT, TANK_NUM) %>%
  reframe(NEP_mean = mean(NEP, na.rm = TRUE))

NEP_means <- NEP_sum_data %>% 
  group_by(TREATMENT) %>% 
  summarize(NEP_total_mean = mean(NEP_mean, na.rm = TRUE),
            NEP_se = sd(NEP_mean, na.rm = TRUE)/sqrt(n()))
NEP_means

## mean NEP per treatment plot #
NEP_mean_plot <- NEP_means %>%
  ggplot(aes(x = TREATMENT, y = NEP_total_mean, color = TREATMENT)) +
  labs(x = "Community Tank", y = expression(bold("Mean Net Ecosystem Production" ~ (mmol ~ C)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble/CCA-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = NEP_sum_data, aes(x = TREATMENT, y = NEP_mean), alpha = 0.7) +
  geom_errorbar(aes(ymin = NEP_total_mean - NEP_se, ymax = NEP_total_mean + NEP_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEP_mean_plot
#ggsave(plot = NEP_mean_plot, filename = here("Output", "NEP_means_plot.png"), width = 9, height = 9)

###################################
# STATS # 

## GROSS PHOTOSYNTHESIS ## 

# how is GP impacted by treatment type # 
GP_treatment_model <- lmer(GP ~ TREATMENT + (1|GENOTYPE), data = final_respo_meta_join)
check_model(GP_treatment_model)
summary(GP_treatment_model)
anova(GP_treatment_model) # the above model shows that the algae and rubble dominated communities had the most significant effect on altering 
# GP of corals in the response tanks (p < 0.05) 
emmeans(GP_treatment_model, pairwise ~ "TREATMENT", adjust="Tukey") # emmeans: estimated marginal means
# summary of average response across levels of TREATMENT to compare means and interpret how each 
# TREATMENT affects GP 
# see most significant p value at contrast between control - algae dom --> which I believe means that the 
# algae dom community is most significantly different from the control? 

# mean pH effect on GP #
GP_pH_model <- lmer(GP ~ pH_mean + (1|GENOTYPE), data = final_respo_meta_join)
check_model(GP_pH_model)
summary(GP_pH_model) # no significant effect of mean pH on coral GP

# range in pH effect on GP # 
GP_pHrange_model <- lmer(GP ~ pH_range + (1|GENOTYPE), data = final_respo_meta_join)
check_model(GP_pHrange_model)
summary(GP_pHrange_model) # no significant effect of pH range on coral GP 

## RESPIRATION ## 

R_treatment_model <- lmer(R ~ TREATMENT + (1| GENOTYPE), data = final_respo_meta_join)
check_model(R_treatment_model)
summary(R_treatment_model) # sig effect of algae dom on coral R (p = 0.035) 
# sig effect of rubble dom on coral R (p = 0.0055) 
# corals in response tanks are most impacted by algae and rubble dominated tanks -- what are the DOC contents of these tanks?
# could the DOC in these communities be altering these respiration rates?
anova(R_treatment_model) # sig effect of community type on R (p < 0.05) 

emmeans(R_treatment_model, pairwise ~ "TREATMENT", adjust="Tukey")
# most significant p values in contrasts of control - rubble, and coral - rubble 

# effect of mean pH on coral R #
R_pHmean_model <- lm(R ~ pH_mean, data = final_respo_meta_join)
check_model(R_pHmean_model)
summary(R_pHmean_model) # sig effect of mean pH on coral R rate (p = 0.0310)

# effect of range in pH on coral R # 
R_pHrange_model <- lm(R ~ pH_range, data = final_respo_meta_join)
check_model(R_pHrange_model)
summary(R_pHrange_model) # sig effect of pH range on coral R (p < 0.05)

## from the above analyses, it looks like respiration rates in the coral response tanks were the most affected by pH
## and community type. need to look further into impact of TA and DOC


## NET PHOTOSYNTHESIS ## 

NP_treatment_model <- lmer(NP ~ TREATMENT + (1|GENOTYPE), data = final_respo_meta_join)
check_model(NP_treatment_model)
summary(NP_treatment_model) # algae dom community had most significant effect on differences in coral NP 
anova(NP_treatment_model) # however, no significant effect of community type? 
emmeans(R_treatment_model, pairwise ~ "TREATMENT", adjust="Tukey")
 
# effect of mean pH on coral response NP # 
NP_pHmean_model <- lm(NP ~ pH_mean, data = final_respo_meta_join)
check_model(NP_pHmean_model)
summary(NP_pHmean_model) # no sig effect of mean pH on coral NP 

# effect of range in pH on coral response NP # 
NP_pHrange_model <- lm(NP ~ pH_range, data = final_respo_meta_join)
check_model(NP_pHrange_model)
summary(NP_pHrange_model) # no sig effect of mean pH on coral NP 


## GP/R Ratio ## 

GPR_treatment_model <- lmer(GP/R ~ TREATMENT + (1|GENOTYPE), data = final_respo_meta_join)
check_model(GPR_treatment_model)
summary(GPR_treatment_model)
anova(GPR_treatment_model) # no sig effect of community type on GP/R ratio

emmeans(R_treatment_model, pairwise ~ "TREATMENT", adjust="Tukey")
# most significant effect in contrasts of control - rubble and coral - rubble 

ggplot(data = final_respo_meta_join, aes(x = TREATMENT, y = pH_mean)) +
  geom_jitter()
anova(lm(pH_range ~ TREATMENT, data = final_respo_meta_join))


## NEP ## 
NEP_treatment_model <- lmer(NEP ~ TREATMENT + (1|GENOTYPE) + (1|TANK_NUM), data = physio_metadata)
check_model(NEP_treatment_model)
summary(NEP_treatment_model)
anova(NEP_treatment_model)

### EFFECTS OF DOC AND CORAL METABOLISM ### 

# read in full DOC data 
DOC_full <- read_csv(here("Data", "DOC", "DOC_full_data.csv"))
# join DOC and final respo data sheets 
final_respo_meta_join$TANK_NUM <- as.character(final_respo_meta_join$TANK_NUM)

final_respo_meta_DOC <- final_respo_meta_join %>% 
  right_join(DOC_full) %>% 
  select(DATETIME, TANK_NUM, GENOTYPE, TREATMENT, GP, R, NP, pH_mean, pH_range, NPOC_mg_L, TN_mg_L, NPOC_uM, TN_uM)



