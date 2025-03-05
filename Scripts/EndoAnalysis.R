library(tidyverse)
library(here)
library(moments)
library(emmeans)
library(agricolae)

here()
endos <- read_csv(here("Data", "Endosymbionts", "MOBEAST_FCM_dataframe.csv"))
surface_area <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA_calculated.csv"))
afdw_data <- read_csv(here("Data", "Data_Raw", "Growth", "coral_mean_biomass_calculated.csv"))

## Endosymbiont counts are measured as event/uL

## Calculate mean endo counts per treatment ## 
endo_data <- endos %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, INCLUDE_VOL_UL, ENDO_COUNT) %>%
  group_by(TREATMENT)

## Calculate endo count per cm2 using blastate volume from "AFDW" data set ##
# combine endo_data with afdw data # 
endo_data$CORAL_NUM <- as.numeric(endo_data$CORAL_NUM)
endo_data2 <- endo_data %>%
  right_join(afdw_data) %>%
  right_join(surface_area) %>%
  select(-c(CORALID, weight1_g, weight2_g, weight_of_wax_g, date)) %>% 
  mutate(endo_per_cm2 = (ENDO_COUNT)*(1/1000)*(BLASTATE_VOL_ML)*(1/SA_cm_2)) %>%
  drop_na()

#write_csv(endo_data2, here("Data", "Endosymbionts", "endo_data_calculated.csv"))

endo_plotdata <- endo_data2 %>%
  group_by(TREATMENT) %>%
  summarise(mean_endos = mean(endo_per_cm2, na.rm = TRUE),
            se_endos = sd(endo_per_cm2, na.rm = TRUE)/sqrt(n()))

endo_plotdata$TREATMENT <- factor(endo_plotdata$TREATMENT, levels = c("Control", "Pre", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))
endo_plot <- endo_plotdata %>%
  ggplot(aes(x = TREATMENT, y = mean_endos, color = TREATMENT)) +
  labs(x = "Treatment", y = "Mean Endosymbiont Event per cm2") +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated", 
                            "Pre" = "Pre")) +
  geom_point() +
  geom_jitter(data = endo_data2, aes(x = TREATMENT, y = endo_per_cm2), alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_endos - se_endos,
                    ymax = mean_endos + se_endos), color = "black", width = 0.1) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan", "Pre" = "purple")) +
  theme_classic()
endo_plot

#ggsave(plot = endo_plot, filename = here("Output", "EndoOutput", "mean_endo_per_treatment.png"), width = 10, height = 6)


