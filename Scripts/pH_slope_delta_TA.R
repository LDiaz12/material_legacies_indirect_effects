library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(lubridate)
library(calecopal)
library(ggridges)

## bring in pH calibration files and raw data files
pHcalib<-read_csv(here("Data","TrisCalSummer2024.csv"))
pHData<-read_csv(here("Data", "CarbonateChemistry.csv"))
TableID<-read_csv(here("Data", "TableID.csv"))

pHSlope<-pHcalib %>%
  nest_by(TrisCalDate)%>%
  mutate(fitpH = list(lm(mVTris~TTris, data = data))) %>% # linear regression of mV and temp of the tris
  summarise(broom::tidy(fitpH)) %>% # make the output tidy
  select(TrisCalDate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  right_join(.,pHData) %>% # join with the pH sample data
  mutate(mVTris = TempInLab*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  drop_na(TempInSitu)%>%
  drop_na(mV) %>%
  mutate(pH = pH(Ex=mV,Etris=mVTris,S=Salinity,T=TempInLab))  # calculate pH of the samples using the pH seacarb function

#Now calculate pH
pHSlope <-pHSlope%>%
  mutate(pH_insitu = pHinsi(pH = pH, ALK = 2200, Tinsi = TempInSitu, Tlab = TempInLab, 
                            S = Salinity,Pt = 0.1, k1k2 = "m10", kf = "dg")) %>%
  select(!pH)%>%
  rename(pH = pH_insitu) %>% # rename it 
  ungroup() %>%
  select(-c(TempInLab, mV, TrisCalDate, TTris, `(Intercept)`, mVTris))

# remove the inflow data and join it with the tanks that had that specific inflow water
InflowData <- pHSlope %>%
  filter(TankID %in% c("Inflow1","Inflow2")) %>%
  select(-c(Flow_Left_30s, Flow_Right_30s, Notes, DO_mg_L, DO_Percent, Salinity, TempInSitu))  %>% ### remove the values that I don't need -- You will eventually need to keep TA which is why I dropped these instead of coding for the ones that I need
  rename(pH_inflow = pH,
         TA_inflow = TA) %>%# rename the pH to show that it is inflow pH
  mutate(InflowTable = ifelse(TankID == "Inflow1",1,2)) %>% # give them inflow numbers to pair easilty with the TankID 
  ungroup()%>%
  select(Date,Time,InflowTable, pH_inflow, TA_inflow) # drop the Tank ID column to be able to join the data correctly by inflow #

SurfaceArea <- 22.5*22.5 # put in the surface area in cm2 for the bottom of the tank here

Data<-pHSlope %>%
  ungroup()%>%
  filter(!TankID %in% c("Inflow1","Inflow2"))%>% # filter out the inflow data now
  mutate(TankID = as.numeric(TankID))%>% # covert to numberics since the inflow data is now dropped
  left_join(TableID) %>%
  left_join(InflowData) %>%# join with the inflowdata for easier calculations of rates
  mutate(DateTime = mdy_hms(paste(Date,Time)),# make a datetime
         pHDiff = pH - pH_inflow,# calculate the difference between the inflow and the pH in each tank 
         totalflow = Flow_Right_30s+Flow_Left_30s,
         residence_time = (1/totalflow)*(10000/60),# convert ml/min to hours by multiplying by the volumner of water in ml and divide by 60
         deltaTA = TA_inflow - TA, # calculate the difference between in and outflow
         NEC = (deltaTA/2)*(1.025)*(10)*(1/residence_time)*(1/SurfaceArea) ### for a real rate should probably normalize the delta TA to the delta control just like in respo
         ) 

### Now Make a plot showing how the Tank pH differed from the inflow pH over time

tank_pH_diffs <- Data %>%
  ggplot(aes(x = DateTime, y = pHDiff, color = Treatment, group = TankID))+
  geom_point()+
  geom_line()
tank_pH_diffs

# Now do the average
avg_pH_treatment_time <- Data %>%
  group_by(Treatment, DateTime)%>%
  summarise(mean_diff = mean(pHDiff, na.rm = TRUE),
            se_diff = sd(pHDiff, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = DateTime, y = mean_diff, color = Treatment))+
  geom_rect(aes(xmin = ymd_hms("2024-06-02 11:30:00"), xmax = ymd_hms("2024-06-02 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-03 06:00:00"), xmax = ymd_hms("2024-06-03 12:30:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-03 12:30:00"), xmax = ymd_hms("2024-06-03 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-03 18:00:00"), xmax = ymd_hms("2024-06-04 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-02 18:00:00"), xmax = ymd_hms("2024-06-03 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+ ## add colors for light and dark times
  geom_rect(aes(xmin = ymd_hms("2024-06-04 06:00:00"), xmax = ymd_hms("2024-06-04 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-04 18:00:00"), xmax = ymd_hms("2024-06-05 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-05 06:00:00"), xmax = ymd_hms("2024-06-05 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-05 18:00:00"), xmax = ymd_hms("2024-06-06 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_hline(yintercept = 0, lty = 2)+ # show where values shifts from positive to negative
  geom_point(size = 1.5)+
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff+se_diff), width = 0.1)+
  geom_line()+
  annotate("text", x = ymd_hms("2024-06-03 08:00:00"), y = 0.1, label = "Overcast", size = 5)+
  labs(x="",
       y = "Change in pH due to community")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
avg_pH_treatment_time

## deltaTA plots
delta_TA <- Data %>%
  ggplot(aes(x = factor(DateTime), y = NEC, color = Treatment))+
  geom_hline(aes(yintercept = 0))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~Treatment)
delta_TA
