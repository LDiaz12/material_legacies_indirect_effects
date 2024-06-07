rm(list=ls())

## process pH

library(tidyverse)
library(seacarb)
library(broom)
library(lubridate)
library(here)
here()
## bring in pH calibration files and raw data files
pHcalib<-read_csv(here("Data","TrisCalSummer2024.csv"))

#pHcalib<-read_csv('~/Desktop/Repositories/material_legacies_indirect_effects/Data/TrisCal05_28.csv')
pHData<-read_csv(here("Data", "CarbonateChemistry.csv"))
#pHData<-read_csv('~/Desktop/Repositories/material_legacies_indirect_effects/Data/Community_Test_Measurements.csv')

## take the mV calibration files by each date and use them to calculate pH

pHSlope<-pHcalib %>%
  nest_by(TrisCalDate)%>%
  mutate(fitpH = list(lm(mVTris~TTris, data = pHcalib))) %>% # linear regression of mV and temp of the tris
  reframe(broom::tidy(fitpH)) %>% # make the output tidy
  select(TrisCalDate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  left_join(.,pHData) %>% # join with the pH sample data
  mutate(mVTris = TempInLab*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  mutate(pH = pH(Ex=mV,Etris=mVTris,S=Salinity,T=TempInLab)) %>% # calculate pH of the samples using the pH seacarb function
  #mutate(pH_insitu = pHinsi(pH = pH, ALK = TA_Raw, Tinsi = TempInSitu, Tlab = Temp, S = Salinity_lab_Silbiger)) %>%
  mutate(pH_insi = pHinsi(pH = pH, Tinsi = TempInSitu, Tlab = TempInLab, S = Salinity, pHscale = "TRUE" )) %>%
  select(Date, Rep, Treatment, TankID, Salinity,pH = pH_insi, TempInSitu, DO, DO_mg_L, Time) %>% ## need to calculate pH insi
  mutate(DateTime = paste(Date, Time),
         DateTime = ymd_hms(DateTime))
  
  
View(pHSlope)

## write the data
# update daily
write_csv(x = pHSlope, file = here("Data", "pH_Slope_06_05.csv"))


pH_plot <- pHSlope %>%
  group_by(Treatment, DateTime, Date) %>%
  summarise(meanpH = mean(pH, na.rm = TRUE), 
            meanDO = mean(DO, na.rm = TRUE), 
            meanDO_mg_L = mean(DO_mg_L, na.rm = TRUE), 
            sepH = sd(pH, na.rm = TRUE)/sqrt(n()), 
            seDO = sd(DO, na.rm = TRUE)/sqrt(n()), 
            seDO_mg_L = sd(DO_mg_L, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = DateTime, y = meanpH, color = Treatment, group = Treatment)) + 
  geom_point() +
  geom_errorbar( aes(ymin = meanpH - sepH, ymax = meanpH + sepH), width = 0.1)+
  geom_line() +
  facet_wrap(~Date, scale = "free")

ggsave(plot = pH_plot, filename = here("Output", "pH_plot.png"), width = 9, height = 6)


ggplot(pHSlope, aes(x = DateTime, y = DO_mg_L, color = Treatment, group = TankID))+
  geom_point() +
  geom_line() +
  #geom_label(aes(label = TankID)) +
  facet_wrap(~Date, scale = "free")
  

