library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(lubridate)
library(ggridges)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(moments)
library(emmeans)
library(agricolae)
library(tidyr)
library(car)


# bring in all growth files 
diaz_std_curve<-read_csv(here("Data_Raw","Growth", "SA", "Laurel_Standard_Curve.csv"))
surface_area<-read_csv(here("Data_Raw","Growth", "SA", "MO24BEAST_SA.csv"))
bw_le<-read_csv(here("Data_Raw","Growth", "PRUS_BW_LE.csv"))

# clean up growth data and join into one dataframe
bw_le_clean <- bw_le %>%
  select(CoralID, Genotype, Treatment, Initial_weight_g, Final_weight_g, Initial_LE_mm, Final_LE_mm)
sa_clean <- surface_area %>%
  select(CoralID, weight1_g, weight2_g, weight_of_wax_g, SA_cm_2)

growth_data <- right_join(bw_le_clean, sa_clean, by='CoralID', relationship = "many-to-many") %>%
  drop_na(Final_weight_g, Initial_weight_g, Final_LE_mm, Initial_LE_mm)

Num_days <- 25

# calculate total growth and LE and add column to dataframe
growth_data <- growth_data %>%
  mutate(Total_Growth = Final_weight_g - Initial_weight_g,
         Total_LE = Final_LE_mm - Initial_LE_mm,
         Growth_Rate = (Total_Growth/(SA_cm_2*Num_days))) # growth rate normalized to g per cm2 per 25 day experimental period



# check assumptions are met
skewness(growth_data$Growth_Rate)
kurtosis(growth_data$Growth_Rate)
hist(growth_data$Growth_Rate)

skewness(growth_data$Total_LE) 
kurtosis(growth_data$Total_LE)
hist(growth_data$Total_LE)

# ANOVAs for BW and LE - include random effect of genotype
growth_data$Genotype <- as.factor(growth_data$Genotype)

bw_geno_model <- lmer(Growth_Rate~Treatment + (1|Genotype) + (1|Genotype:Treatment), data=growth_data)
qqp(residuals(bw_geno_model), "norm")
summary(bw_geno_model)
anova(bw_geno_model)


le_geno_model <- lmer(Total_LE ~ Treatment + (1|Genotype) + (1|Genotype:Treatment), data=growth_data)
qqp(residuals(le_geno_model), "norm")
#le_log <- log(growth_data$LE_Normalized)
summary(le_geno_model)
anova(le_geno_model)


# simple box plot visualization
#bw_plot <- ggplot(growth_data, aes(x=Treatment, y=Growth_Normalized, fill=Treatment)) +
  #geom_boxplot()
#bw_plot
bw_gen_plot <- ggplot(growth_data, aes(x=Treatment, y=Growth_Rate, fill=Treatment)) + 
  geom_boxplot() + 
  labs(x="Treatment", y="Growth (g/d-1cm-2)") + 
  scale_fill_manual(values=c("Control"="green","Algae-Dom"="red","Coral-Dom" = "blue", "Rubble-Dom" = "yellow"), name="Treatment") +
  theme_bw(base_size=14) +
  geom_point(position = position_dodge(width=0.4))
bw_gen_plot

#le_plot <- ggplot(growth_data, aes(x=Treatment, y=LE_Normalized, fill=Treatment)) +
  #geom_boxplot()
#le_plot

le_gen_plot <- ggplot(growth_data, aes(x=Treatment, y=Total_LE, fill=Treatment)) + 
  geom_boxplot() + 
  labs(x="Treatment", y="Linear Extension (mm/cm2)") + 
  scale_fill_manual(values=c("Control"="green","Algae-Dom"="red","Coral-Dom" = "blue", "Rubble-Dom" = "yellow"), name="Treatment") +
  theme_bw(base_size=14) +
  geom_point(position = position_dodge(width=0.4))
le_gen_plot
