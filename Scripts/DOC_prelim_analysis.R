library(tidyverse)
library(here)
library(ggplot2)
library(moments)
library(emmeans)
library(agricolae)
library(tidyr)
library(dplyr)

here()
DOC_1 <- read_csv(here("Data", "DOC", "DOC_JUN2_JUN6.csv"))
DOC_2 <- read_csv(here("Data", "DOC", "DOC_JUN6_JUN26.csv"))

## DOC1 plot NPOC mg/l ##
DOC_plot1 <- DOC_1 %>%
  ggplot(aes(x = DATE, y = NPOC_mg_L, color = TREATMENT)) + 
  facet_wrap(~TIME) +
  geom_point()
DOC_plot1

## DOC1 plot TN mg/l ## 
DOC_plot1a <- DOC_1 %>%
  ggplot(aes(x = DATE, y = TN_mg_L, color = TREATMENT)) + 
  facet_wrap(~TIME) +
  geom_point()
DOC_plot1a
## DOC1 plot NPOC_uM ##
DOC_plot1b <- DOC_1 %>%
  ggplot(aes(x = DATE, y = NPOC_uM, color = TREATMENT)) + 
  facet_wrap(~TIME) +
  geom_point()
DOC_plot1b
## DOC1 plot TN_uM ## 
DOC_plot1c <- DOC_1 %>%
  ggplot(aes(x = DATE, y = TN_uM, color = TREATMENT)) + 
  facet_wrap(~TIME) +
  geom_point()
DOC_plot1c


## DOC2 plot NPOC mg/L ##
DOC_plot2 <- DOC_2 %>%
  ggplot(aes(x = DATE, y = NPOC_mg_L, color = TREATMENT)) + 
  facet_wrap(~TIME) +
  geom_point() 
DOC_plot2
## DOC2 plot TN mg/l ## 
DOC_plot2a <- DOC_2 %>%
  ggplot(aes(x = DATE, y = TN_mg_L, color = TREATMENT)) + 
  facet_wrap(~TIME) +
  geom_point()
DOC_plot2a
## DOC2 plot NPOC_uM ##
DOC_plot2b <- DOC_2 %>%
  ggplot(aes(x = DATE, y = NPOC_uM, color = TREATMENT)) + 
  facet_wrap(~TIME) +
  geom_point() 
DOC_plot2b
## DOC2 plot TN_uM ## 
DOC_plot2c <- DOC_2 %>%
  ggplot(aes(x = DATE, y = TN_uM, color = TREATMENT)) + 
  facet_wrap(~TIME) +
  geom_point()
DOC_plot2c

## save all plots 
output_folder <- here("Output","DOCOutput")
ggsave(filename = file.path(output_folder, "DOC_plot1_NPOC_mg_L.png"), plot = DOC_plot1, width = 8, height = 6)
ggsave(filename = file.path(output_folder, "DOC_plot1a_TN_mg_L.png"), plot = DOC_plot1a, width = 8, height = 6)
ggsave(filename = file.path(output_folder, "DOC_plot1b_NPOC_uM.png"), plot = DOC_plot1b, width = 8, height = 6)
ggsave(filename = file.path(output_folder, "DOC_plot1c_TN_uM.png"), plot = DOC_plot1c, width = 8, height = 6)

ggsave(filename = file.path(output_folder, "DOC_plot2_NPOC_mg_L.png"), plot = DOC_plot2, width = 10, height = 6)
ggsave(filename = file.path(output_folder, "DOC_plot2a_TN_mg_L.png"), plot = DOC_plot2a, width = 10, height = 6)
ggsave(filename = file.path(output_folder, "DOC_plot2b_NPOC_uM.png"), plot = DOC_plot2b, width = 10, height = 6)
ggsave(filename = file.path(output_folder, "DOC_plot2c_TN_uM.png"), plot = DOC_plot2c, width = 10, height = 6)








