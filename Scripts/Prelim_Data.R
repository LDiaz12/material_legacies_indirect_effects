library(tidyverse)
library(broom)
library(seacarb)

pHcalib<-read_csv('~/Desktop/Prelim_Jan_Data/Tris_Cal_Log_Jan_2024.csv')
pHData<-read_csv('~/Desktop/Prelim_Jan_Data.csv')

## take the mV calibration files by each date and use them to calculate pH
pHSlope<-pHcalib %>%
  group_by(TrisCalDate)%>%
  do(fitpH = lm(mVTris~TTris, data = pHcalib))%>% # linear regression of mV and temp of the tris
  tidy(fitpH) %>% # make the output tidy
  select(SamplingDate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  left_join(.,pHData) %>% # join with the pH sample data
  mutate(mVTris = TempInLab*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  mutate(pH = pH(Ex=mV,Etris=mVTris,S=Salinity,T=TempInLab)) %>% # calculate pH of the samples using the pH seacarb function
  #mutate(pH_insitu = pHinsi(pH = pH, ALK = TA_Raw, Tinsi = TempInSitu, Tlab = Temp, S = Salinity_lab_Silbiger)) %>%
  select(TrisCalDate, SampleID,Salinity,pH, TempInSitu, TempInLab) ## need to calculate pH insi then it is done

pHSlope <-pHSlope%>%
  mutate(pH_insitu = pHinsi(pH = pH, Tinsi = TempInSitu, Tlab = TempInLab, 
                            S = Salinity,Pt = 0, k1k2 = "m10", kf = "dg")) %>%
  select(!pH) %>% # I only need the in situ pH calculation so remove this
  rename(pH = pH_insitu) %>% # rename it 
  ungroup()


# or
pHSlope<-pHcalib %>%
  nest_by(TrisCalDate)%>%
  mutate(fitpH = list(lm(mVTris~TTris, data = pHcalib))) %>% # linear regression of mV and temp of the tris
  reframe(broom::tidy(fitpH)) %>% # make the output tidy
  select(TrisCalDate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  cross_join(.,pHData) %>% # join with the pH sample data
  mutate(mVTris = TempInSitu*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  mutate(pH = pH(Ex=mV,Etris=mVTris,S=Salinity,T=TempInSitu)) %>% # calculate pH of the samples using the pH seacarb function
  #mutate(pH_insitu = pHinsi(pH = pH, ALK = TA_Raw, Tinsi = TempInSitu, Tlab = Temp, S = Salinity_lab_Silbiger)) %>%
  select(SamplingDate, UniqueID,Salinity,pH, TempInSitu, DO, TempInLab) ## need to calculate pH insi then it is done
pHSlope <-pHSlope%>%
  mutate(pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
                            S = Salinity,Pt = 0, k1k2 = "m10", kf = "dg")) %>%
  select(!pH) %>% # I only need the in situ pH calculation so remove this
  rename(pH = pH_insitu) %>% # rename it 
  ungroup()


View(pHSlope)

## write the data
write.csv(x = pHSlope, file = '~/Desktop/Prelim_Jan_Data.csv')
