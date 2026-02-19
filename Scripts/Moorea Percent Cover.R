library(ggplot2)
library(dplyr)
library(grid)
library(tidyverse)
library(here)

# using LTER1 backreef data
backreef_mean_cover <- read_csv(here("Data", "PercentCoverData", "backreef_benthic_comp_mean_perc.csv"))

# create percent cover plot by benthic category over time
year_fact = cut(backreef_mean_cover$year, 3, labels = c("2006", "2010", "2019"))
table(year_fact)
is.factor(year_fact)

backreef_mean_cover$benthic_group <- factor(backreef_mean_cover$benthic_group, 
                                            levels = c("coral", "ctb", "CCA", "macroalgae", "sand"),
                                            labels = c("Coral", "Rubble", "Crustose Coralline Algae", "Macroalgae", "Sand"))


br_benthic_cover_plot <- ggplot(data=backreef_mean_cover, aes(x=year_fact, y=mean_percent_cover*100, fill=benthic_group)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Year", y = "Mean Percent Cover", fill = "Benthic Category") +
  scale_fill_manual(values = c("Coral" = "coral", "Rubble" = "lightgray", "Crustose Coralline Algae" = "orchid" ,"Macroalgae" = "darkgreen", "Sand" = "tan"),
                    labels = c("Coral", "Rubble", "Crustose Coralline Algae", "Macroalgae", "Sand")) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = "bold")) +
  geom_text(aes(label = paste0(round(mean_percent_cover, 0), "%")),
            position = position_fill(vjust = 0.5), # centers text within each bar segment
            color = "black", size = 4)
br_benthic_cover_plot
#ggsave(plot = br_benthic_cover_plot, filename = here("Output", "br_benthic_cover_plot.pdf"), width = 6, height = 6)
