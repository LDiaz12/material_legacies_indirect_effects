library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(lubridate)
library(ggridges)

## bring in pH calibration files and raw data files
pHcalib<-read_csv(here("Data","TrisCalSummer2024.csv"))
pHData<-read_csv(here("Data", "CarbonateChemistry.csv"))
TableID<-read_csv(here("Data", "TableID.csv"))

## update daily!
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
  select(-c(Flow_Left_30s, Flow_Right_30s, Notes, DO_mg_L, Salinity, TempInSitu))  %>% ### remove the values that I don't need -- You will eventually need to keep TA which is why I dropped these instead of coding for the ones that I need
  rename(pH_inflow = pH,
         TA_inflow = TA) %>%# rename the pH to show that it is inflow pH
  mutate(InflowTable = ifelse(TankID == "Inflow1",1,2)) %>% # give them inflow numbers to pair easily with the TankID 
  ungroup()%>%
  select(Date,Time,InflowTable, pH_inflow, TA_inflow) # drop the Tank ID column to be able to join the data correctly by inflow #

SurfaceArea <- 22.5*22.5 # put in the surface area in cm2 for the bottom of the tank here

Data<-pHSlope %>%
  ungroup()%>%
  filter(!TankID %in% c("Inflow1","Inflow2"))%>% # filter out the inflow data now
  mutate(TankID = as.numeric(TankID))%>% # convert to numeric since the inflow data is now dropped
  left_join(TableID) %>%
  left_join(InflowData) %>% # join with the inflow data for easier calculations of rates
  mutate(DateTime = ymd_hms(paste(Date,Time)), # make a datetime
         pHDiff = pH - pH_inflow, # calculate the difference between the inflow and the pH in each tank 
         totalflow = Flow_Right_30s+Flow_Left_30s,
         residence_time = (1/totalflow)*(10000/60),# convert ml/min to hours by multiplying by the volumne of water in ml and divide by 60
         deltaTA = TA_inflow - TA, # calculate the difference between in and outflow
         NEC = (deltaTA/2)*(1.025)*(10)*(1/residence_time)*(1/SurfaceArea) ### for a real rate should probably normalize the delta TA to the delta control just like in respo
         )

### Now Make a plot showing how the Tank pH differed from the inflow pH over time

tank_pH_diffs <- Data %>%
  ggplot(aes(x = DateTime, y = pHDiff, color = Treatment, group = TankID))+
  geom_point()+
  geom_line()
tank_pH_diffs +
  scale_color_hue(labels = c("Algae-Dominated", "Control", "Rubble-Dominated", "Coral-Dominated"))
#ggsave(plot = tank_pH_diffs, filename = here("Output", "tank_pH_diffs.png"), width = 11, height = 9)

# Now do the average
avg_pH_treatment_time <- Data %>%
  group_by(Treatment, DateTime)%>%
  summarise(mean_diff = mean(pHDiff, na.rm = TRUE),
            se_diff = sd(pHDiff, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = DateTime, y = mean_diff, color = Treatment))+
  geom_rect(aes(xmin = ymd_hms("2024-06-02 06:00:00"), xmax = ymd_hms("2024-06-02 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-02 18:00:00"), xmax = ymd_hms("2024-06-03 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-03 06:00:00"), xmax = ymd_hms("2024-06-03 12:30:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-03 12:30:00"), xmax = ymd_hms("2024-06-03 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-03 18:00:00"), xmax = ymd_hms("2024-06-04 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+ 
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
            fill = "lightgrey", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-06 06:00:00"), xmax = ymd_hms("2024-06-06 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-06 18:00:00"), xmax = ymd_hms("2024-06-07 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-07 06:00:00"), xmax = ymd_hms("2024-06-07 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-07 18:00:00"), xmax = ymd_hms("2024-06-08 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-08 06:00:00"), xmax = ymd_hms("2024-06-08 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-08 18:00:00"), xmax = ymd_hms("2024-06-09 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-09 06:00:00"), xmax = ymd_hms("2024-06-09 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-09 18:00:00"), xmax = ymd_hms("2024-06-10 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-10 06:00:00"), xmax = ymd_hms("2024-06-10 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-10 18:00:00"), xmax = ymd_hms("2024-06-11 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-11 06:00:00"), xmax = ymd_hms("2024-06-11 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-11 18:00:00"), xmax = ymd_hms("2024-06-12 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-12 06:00:00"), xmax = ymd_hms("2024-06-12 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-12 18:00:00"), xmax = ymd_hms("2024-06-13 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-13 06:00:00"), xmax = ymd_hms("2024-06-13 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-13 18:00:00"), xmax = ymd_hms("2024-06-14 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-14 06:00:00"), xmax = ymd_hms("2024-06-14 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-14 18:00:00"), xmax = ymd_hms("2024-06-15 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-15 06:00:00"), xmax = ymd_hms("2024-06-15 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-15 18:00:00"), xmax = ymd_hms("2024-06-16 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-16 06:00:00"), xmax = ymd_hms("2024-06-16 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-16 18:00:00"), xmax = ymd_hms("2024-06-17 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-17 06:00:00"), xmax = ymd_hms("2024-06-17 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-17 18:00:00"), xmax = ymd_hms("2024-06-18 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-18 06:00:00"), xmax = ymd_hms("2024-06-18 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-18 18:00:00"), xmax = ymd_hms("2024-06-19 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-19 06:00:00"), xmax = ymd_hms("2024-06-19 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-19 18:00:00"), xmax = ymd_hms("2024-06-20 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-20 06:00:00"), xmax = ymd_hms("2024-06-21 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-21 18:00:00"), xmax = ymd_hms("2024-06-22 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-22 06:00:00"), xmax = ymd_hms("2024-06-22 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-22 18:00:00"), xmax = ymd_hms("2024-06-23 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-23 06:00:00"), xmax = ymd_hms("2024-06-23 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-23 18:00:00"), xmax = ymd_hms("2024-06-24 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-24 06:00:00"), xmax = ymd_hms("2024-06-24 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-24 18:00:00"), xmax = ymd_hms("2024-06-25 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-25 06:00:00"), xmax = ymd_hms("2024-06-25 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-25 18:00:00"), xmax = ymd_hms("2024-06-26 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-26 06:00:00"), xmax = ymd_hms("2024-06-26 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-26 18:00:00"), xmax = ymd_hms("2024-06-27 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_hline(yintercept = 0, lty = 2)+ # show where values shifts from positive to negative
  geom_point(size = 1.5)+
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff+se_diff), width = 0.1)+
  geom_line()+
  annotate("text", x = ymd_hms("2024-06-03 08:00:00"), y = 0.1, label = "Overcast", size = 5)+
  annotate("text", x = ymd_hms("2024-06-06 13:00:00"), y = 0.1, label = "Overcast \n Rain", size = 5)+
  labs(x="",
       y = "Change in pH due to community")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
avg_pH_treatment_time +
  scale_color_hue(labels = c("Algae-Dominated", "Control", "Rubble-Dominated", "Coral-Dominated"))
ggsave(plot = avg_pH_treatment_time, filename = here("Output", "avg_pH_treatment_time.png"), width = 15, height = 10)

## calculating total flow and residence time. add pHdiff and deltaTA
Data<-pHSlope %>%
  ungroup()%>%
  filter(!TankID %in% c("Inflow1","Inflow2"))%>% # filter out the inflow data now
  mutate(TankID = as.numeric(TankID))%>% # convert to numeric since the inflow data is now dropped
  left_join(TableID) %>%
  left_join(InflowData) %>% # join with the inflow data for easier calculations of rates
  mutate(DateTime = ymd_hms(paste(Date,Time)), # make a datetime
         pHDiff = pH - pH_inflow, # calculate the difference between the inflow and the pH in each tank 
         totalflow = Flow_Right_30s+Flow_Left_30s,
         residence_time = (1/totalflow)*(11356.2/60), # convert ml/min to hours by multiplying by the volume of water in ml and divide by 60
         deltaTA = TA_inflow - TA) # calculate the difference between in and outflow

# pull out the delta TA from the controls and take average by inflow and normalize the NEC rates to it-- this accounts for changes in TA due to background water
control_deltaTA<-Data %>%
  filter(Treatment == "Control") %>%
  select(DateTime, deltaTA,InflowTable ) %>%
  group_by(DateTime, InflowTable)%>%
  summarise(deltaTA_blank = mean(deltaTA, na.rm = TRUE))


Data<-Data %>%
  left_join(control_deltaTA) %>%
  mutate(NEC = ((deltaTA-deltaTA_blank)/2)*(1.025)*(10)*(1/residence_time)*(1/SurfaceArea) ### for a real rate should probably normalize the delta TA to the delta control just like in respo
  )

delta_TAs <- Data %>%
  ggplot(aes(x = DateTime, y = deltaTA, color = Treatment, na.rm = TRUE))+
  geom_point()+
  geom_line()
delta_TAs +
  scale_color_hue(labels = c("Algae-Dominated", "Control", "Rubble-Dominated", "Coral-Dominated"))
ggsave(plot = delta_TAs, filename = here("Output", "delta_TAs.png"), width = 10, height = 9)


# light vs diff pH plot
light_pH <- Data %>%
  ggplot(aes(x = Light_nm, y = pHDiff, color = Treatment))+
  geom_point()+
  geom_line()
light_pH +
  scale_color_hue(labels = c("Algae-Dominated", "Control", "Rubble-Dominated", "Coral-Dominated"))
ggsave(plot = light_pH, filename = here("Output", "light_pH.png"), width = 15, height = 9)
