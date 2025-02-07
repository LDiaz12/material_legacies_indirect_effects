library(tidyverse)
library(here)
library(moments)
library(emmeans)
library(agricolae)

here()
endos <- read_csv(here("Data", "Endosymbionts", "MOBEAST_FCM_dataframe.csv"))
surface_area <- read_csv(here("Data", "Data_Raw", "Growth", "SA", "MO24BEAST_SA_calculated.csv"))
afdw_data <- read_csv(here("Data", "Data_Raw", "Growth", "MO24BEAST_AFDW.csv"))

## Endosymbiont counts are measured as event/uL

## Calculate mean endo counts per treatment ## 
endo_data <- endos %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, INCLUDE_VOL_UL, ENDO_COUNT) %>%
  group_by(TREATMENT) %>%
  #mutate(mean_endos = mean(EndosymbiontCount, na.rm = TRUE)) %>% ## mean endo counts per treatment
  mutate()
  drop_na()

## Calculate endo count per cm2 using blastate volume from "AFDW" data set ##



  
endo_data$Treatment <- factor(endo_data$Treatment, levels = c("BLANK", "Control", "Pre", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))
endo_plot <- endo_data %>%
  group_by(Treatment) %>%
  summarise(mean_endos = mean(EndosymbiontCount, na.rm = TRUE),
            se_endos = sd(EndosymbiontCount, na.rm = TRUE)/sqrt(n())) %>%
  ggplot(aes(x = Treatment, y = mean_endos, color = Treatment, na.rm = TRUE)) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated", 
                            "Pre" = "Pre", "BLANK" = "Blank")) +
  scale_fill_identity() +
  geom_point(size = 2.5) +
  geom_jitter(data = endo_data, aes(x = Treatment, y = EndosymbiontCount), alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_endos - se_endos,
                    ymax = mean_endos + se_endos), color = "black", width = 0.1) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan", "Pre" = "purple", "Blank" = "black")) +
  theme_classic()
endo_plot

output_folder <- here("Output","EndoOutput")
ggsave(filename = file.path(output_folder, "endo_plot.png"), plot = endo_plot, width = 10, height = 6)


