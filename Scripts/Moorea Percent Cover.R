library(ggplot2)
library(dplyr)
library(grid)
library(tidyverse)
library(here)
here()
# determining most dominant coral spp on the backreef of LTER1
coral_cover_br <- read_csv(here("Data", "PercentCoverData", "corals.csv"))
coral_cover_br_plot <- ggplot(data = coral_cover_br, aes(fill=name, y=mean_value, x=YEAR)) +
  geom_bar(position="fill", stat="identity") +
  labs(x="Year", y="Mean Percent Cover", title = "Backreef Coral Cover", fill = "Coral Species")
coral_cover_br_plot
# plot shows clear dominance of porites spp across disturbance time points

# using LTER1 backreef data
backreef_mean_cover <- read_csv(here("Data", "PercentCoverData", "backreef_benthic_comp_mean_perc.csv"))

# create percent cover plot by benthic category over time
year_fact = cut(backreef_mean_cover$year, 3, labels = c("2006", "2010", "2019"))
table(year_fact)
is.factor(year_fact)

backreef_mean_cover$benthic_group <- factor(backreef_mean_cover$benthic_group, 
                                            levels = c("coral", "ctb", "macroalgae", "sand"),
                                            labels = c("Coral", "Rubble", "Macroalgae", "Sand"))


br_benthic_cover_plot <- ggplot(data=backreef_mean_cover, aes(fill=benthic_group, x=year_fact, y=mean_percent_cover)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Year", y = "Mean Percent Cover", title = "Mo'orea Backreef Benthic Percent Cover", fill = "Benthic Category") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = c("Coral" = "coral", "Rubble" = "lightgray", "Macroalgae" = "darkgreen", "Sand" = "tan"),
                    labels = c("Coral", "Rubble", "Macroalgae", "Sand")) +
  geom_text(aes(label = paste0(round(mean_percent_cover, 0), "%")),
            position = position_fill(vjust = 0.5), # centers text within each bar segment
            color = "black", size = 4)
br_benthic_cover_plot
ggsave(plot = br_benthic_cover_plot, filename = here("Output", "br_benthic_cover_plot.png"), width = 6, height = 6)

##pre cover plot##
backreef_2006 <- read_csv(here("Data", "PercentCoverData", "backreef_2006_cover.csv"))
backreef_2006_plot <- ggplot(data=backreef_2006, aes(fill=benthic_category, x=benthic_category, y=percent_cover)) + 
  stat_summary(fun.data=mean_sdl, geom="bar") +
  labs(x="Benthic Category", y="Percent Cover", 
       title = "Pre-Disturbance Percent Cover\nLTER Backreef: Coral Dominated", fill = "Benthic Category") +
  scale_color_fermenter() +
  theme_bw()+
  theme(legend.text = element_text(size=8)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
backreef_2006_plot + theme(axis.text.x = element_text(angle = 40))

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
