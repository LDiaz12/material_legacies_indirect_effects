library(tidyverse)
library(here)
library(moments)
library(emmeans)
library(agricolae)
library(ggrepel)
library(lubridate)

here()
DOC_1 <- read_csv(here("Data", "DOC", "DOC_JUN2_JUN6.csv"))
DOC_2 <- read_csv(here("Data", "DOC", "DOC_JUN6_JUN26.csv"))

DOC_full <- bind_rows(DOC_1, DOC_2) %>%
  mutate(DATE = mdy(DATE))


## DOC1 plot NPOC mg/l ##
DOC_plot1 <- DOC_full %>%
  ggplot(aes(x = TREATMENT, y = NPOC_mg_L, color = as.factor(TIME))) + 
  geom_point() +
  geom_text_repel(aes(label = TANK_ID)) +
  facet_wrap(~DATE) 
DOC_plot1

DOC_plot2 <- DOC_full %>%
  ggplot(aes(x = TREATMENT, y = TN_mg_L, color = as.factor(TIME))) + 
  geom_point() +
  geom_text_repel(aes(label = TANK_ID)) +
  facet_wrap(~DATE) 
DOC_plot2

DOC_plot3 <- DOC_full %>%
  ggplot(aes(x = TREATMENT, y = NPOC_uM, color = as.factor(TIME))) + 
  geom_point() +
  geom_text_repel(aes(label = TANK_ID)) +
  facet_wrap(~DATE) 
DOC_plot3

DOC_plot4 <- DOC_full %>%
  ggplot(aes(x = TREATMENT, y = TN_uM, color = as.factor(TIME))) + 
  geom_point() +
  geom_text_repel(aes(label = TANK_ID)) +
  facet_wrap(~DATE) 
DOC_plot4

## save all plots 
output_folder <- here("Output","DOCOutput")
ggsave(filename = file.path(output_folder, "DOC_NPOC_mg_L.png"), plot = DOC_plot1, width = 12, height = 10)
ggsave(filename = file.path(output_folder, "DOC_TN_mg_L.png"), plot = DOC_plot2, width = 12, height = 10)
ggsave(filename = file.path(output_folder, "DOC_NPOC_uM.png"), plot = DOC_plot3, width = 12, height = 10)
ggsave(filename = file.path(output_folder, "DOC_TN_uM.png"), plot = DOC_plot4, width = 12, height = 10)



