moorea_2005 <- read.csv('~/Desktop/moorea 2005 outer 10.csv')
library(tidyverse)

plot1 <- ggplot(data=moorea_2005, mapping=aes(x=Taxonomy_Substrate_or_Functional_Group, y=Percent_Cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") 
plot1 + theme(axis.text.x = element_text(angle = 90))

moorea_COTS <- read.csv('~/Desktop/moorea 2010 COTS outer 10.csv')
plot2 <- ggplot(data=moorea_COTS, mapping=aes(x=Taxonomy_Substrate_or_Functional_Group, y=Percent_Cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") 
plot2 + theme(axis.text.x = element_text(angle = 90))

moorea_bleaching <- read.csv('~/Desktop/moorea 2015 bleaching outer 10.csv')
plot3 <- ggplot(data=moorea_bleaching, mapping=aes(x=Taxonomy_Substrate_or_Functional_Group, y=Percent_Cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") 
plot3 + theme(axis.text.x = element_text(angle = 90))

#'present' is using 2021 data which is the most recent
moorea_present <- read.csv('~/Desktop/moorea present percent cover.csv')
plot4 <- ggplot(data=moorea_present, mapping=aes(x=Taxonomy_Substrate_or_Functional_Group, y=Percent_Cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") 
plot4 + theme(axis.text.x = element_text(angle = 90))
