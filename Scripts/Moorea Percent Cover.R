backreef_cover <- read.csv("~/Desktop/Master's Thesis/Moorea Percent Cover/Backreef_Percent_Cover.csv")
library(ggplot2)
library(dplyr)

cbpalette <- c("#CC79A7","#999999","#009E73","#E69F00", "#F0E442")

library(grid)
percent_plot <- ggplot(data=backreef_cover, aes(fill=benthic_category, y=percent_cover, x=year)) +
  geom_bar(position="fill", stat="identity") +
  labs(x="\nYear", y="Percent Cover", title = "Back reef percent cover", 
       fill = "Benthic Category") +
  scale_fill_manual(values = cbpalette, name = "Benthic Category", labels = 
                      c("Coral", "Rubble", 
                        "Macroalgae", "Other", "Sand")) +
  theme_bw()
percent_plot

##pre cover plot##
backreef_2005 <- read.csv("~/Desktop/backreef_2005_cover.csv")
backreef_2005_plot <- ggplot(data=backreef_2005, aes(fill=benthic_category, x=benthic_category, y=percent_cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") +
  labs(x="Benthic Category", y="Percent Cover", 
       title = "Pre-Disturbance Percent Cover\nLTER Backreef", fill = "Benthic Category") +
  scale_color_fermenter() +
  theme_bw()+
  theme(legend.text = element_text(size=8)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
backreef_2005_plot + theme(axis.text.x = element_text(angle = 90))
## COTS cover plot##
backreef_2010 <- read.csv("~/Desktop/backreef_2010_cover.csv")
backreef_2010_plot <- ggplot(data=backreef_2010, aes(fill=benthic_category, x=benthic_category, y=percent_cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") +
  labs(x="Benthic Category", y="Percent Cover", 
       title = "Post-COTS Percent Cover\nLTER Backreef", fill = "Benthic Category") +
  scale_color_fermenter() +
  theme_bw()+
  theme(legend.text = element_text(size=8)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
backreef_2010_plot + theme(axis.text.x = element_text(angle = 90))
## bleaching cover plot##
backreef_2019 <- read.csv("~/Desktop/backreef_2019_cover.csv")
backreef_2019_plot <- ggplot(data=backreef_2019, aes(fill=benthic_category, x=benthic_category, y=percent_cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") +
  labs(x="Benthic Category", y="Percent Cover", 
       title = "Post-Bleaching Percent Cover\nLTER Backreef", fill = "Benthic Category") +
  scale_color_fermenter() +
  theme_bw()+
  theme(legend.text = element_text(size=8)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
backreef_2019_plot + theme(axis.text.x = element_text(angle = 90))


backreef_algae <- read.csv("~/Desktop/benthic_algal_cover.csv")
backreef_algae_plot <- ggplot(data=backreef_algae, aes(fill=Taxonomy_Substrate_Functional_Group, x=Taxonomy_Substrate_Functional_Group, y=Percent_Cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") +
  labs(x="Benthic Taxonomy", y="Percent Cover", 
       title = "Benthic Algae Percent Cover\nLTER Backreef", fill = "Benthic Algae") +
  scale_color_fermenter() +
  facet_wrap(~Year)+
  theme_bw()+
  theme(legend.text = element_text(size=8)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
backreef_algae_plot + theme(axis.text.x = element_text(angle = 90))


backreef_algae %>%
  group_by(Site, Year, Taxonomy_Substrate_Functional_Group)%>%
  summarise(mean_perc = mean(Percent_Cover, na.rm = TRUE))
