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
library(vegan)

MOBEAST_metadata <- read_csv(here("Data", "MO24BEAST_Metadata_FULL.csv"))
MOBEAST_metadata <- MOBEAST_metadata %>%
  select(c("endo_per_cm2", "chla.ug.cm2", "mean_tissue_biomass", "R", "NP", "GP", 
           "TA_rangemean", "TA_mean", ""))


