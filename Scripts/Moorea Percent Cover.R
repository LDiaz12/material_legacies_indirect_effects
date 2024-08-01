library(ggplot2)
library(dplyr)
library(grid)
library(tidyverse)
library(here)

# determining most dominant coral spp on the backreef of LTER1
coral_cover_br <- read_csv(here("Data", "corals.csv"))
coral_cover_br_plot <- ggplot(data = coral_cover_br, aes(fill=name, y=mean_value, x=YEAR)) +
  geom_bar(position="fill", stat="identity") +
  labs(x="Year", y="Mean Percent Cover", title = "Backreef Coral Cover", fill = "Coral Species")
coral_cover_br_plot
# plot shows clear dominance of porites spp across disturbance time points

# using LTER1 backreef data
backreef_mean_cover <- read_csv(here("Data", "backreef_benthic_comp_mean_perc.csv"))

# create percent cover plot by benthic category over time
year_fact = cut(backreef_mean_cover$year, 3, labels = c("2006", "2010", "2019"))
table(year_fact)
is.factor(year_fact)

br_benthic_cover_plot <- ggplot(data=backreef_mean_cover, aes(fill=benthic_group, x=year_fact, y=mean_percent_cover)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Year", y = "Mean Percent Cover", title = "LTER 1 Backreef Benthic Percent Cover", fill = "Benthic Category") +
  scale_color_fermenter()+
  theme_bw()
br_benthic_cover_plot + 
  scale_fill_discrete(labels = c("Coral", "CTB (CCA, Turf Algae, \n Bare Space)", "Macroalgae", "Sand"))


#br_benthic_percent_plot <- backreef_cover %>%
  #group_by(benthic_category, year) %>%
  #summarise(mean_perc = mean(percent_cover)) %>%
  #ggplot(aes(x=year, y=mean_perc, color=benthic_category)) +
  #geom_bar(position="fill", stat="identity") +
  #labs(x="Year", y="Percent Cover", title = "Backreef Percent Cover", fill = "Benthic Category") +
  #theme_classic()
#br_benthic_percent_plot

##pre cover plot##
backreef_2006 <- read_csv(here("Data", "backreef_2006_cover.csv"))
backreef_2006_plot <- ggplot(data=backreef_2006, aes(fill=benthic_category, x=benthic_category, y=percent_cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") +
  labs(x="Benthic Category", y="Percent Cover", 
       title = "Pre-Disturbance Percent Cover\nLTER Backreef: Coral Dominated", fill = "Benthic Category") +
  scale_color_fermenter() +
  theme_bw()+
  theme(legend.text = element_text(size=8)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
backreef_2006_plot + theme(axis.text.x = element_text(angle = 90))
## COTS cover plot##
backreef_2011 <- read_csv(here("Data", "backreef_2011_cover.csv"))
backreef_2011_plot <- ggplot(data=backreef_2011, aes(fill=benthic_category, x=benthic_category, y=percent_cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") +
  labs(x="Benthic Category", y="Percent Cover", 
       title = "Post-COTS Percent Cover\nLTER Backreef", fill = "Benthic Category") +
  scale_color_fermenter() +
  theme_bw()+
  theme(legend.text = element_text(size=8)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
backreef_2011_plot + theme(axis.text.x = element_text(angle = 90))
## bleaching cover plot##
backreef_2019 <- read_csv(here("Data", "backreef_2019_cover.csv"))
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

## backreef algal community percent cover ##
backreef_algae <- read_csv(here("Data", "benthic_algal_cover.csv"))
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


backreef_algae_by_year <- backreef_algae %>%
  group_by(Site, Year, Taxonomy_Substrate_Functional_Group)%>%
  summarise(mean_perc = mean(Percent_Cover, na.rm = TRUE))

backreef_algae_by_year
