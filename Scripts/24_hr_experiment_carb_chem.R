library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(lubridate)
library(ggridges)
library(ggplot2)
library(moments)
library(emmeans)
library(agricolae)
library(tidyr)

pHcalib<-read_csv(here("Data","Chemistry", "TrisCalSummer2024.csv"))
pHData<-read_csv(here("Data", "Chemistry", "24_hr_carb_chem.csv"))
TableID<-read_csv(here("Data", "TableID.csv"))

## calculate pH slope using mV ##
pHSlope<-pHcalib %>%
  nest_by(TrisCalDate)%>%
  mutate(fitpH = list(lm(mVTris~TTris, data = data))) %>% # linear regression of mV and temp of the tris
  summarise(broom::tidy(fitpH)) %>% # make the output tidy
  select(TrisCalDate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  right_join(.,pHData) %>% # join with the pH sample data
  mutate(mVTris = TEMPINLAB*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  drop_na(TEMPINSITU)%>%
  drop_na(mV) %>%
  mutate(pH = seacarb::pH(Ex=mV,Etris=mVTris,S=SALINITY,T=TEMPINLAB))  # calculate pH of the samples using the pH seacarb function

## calculate pH insitu ## 
pHSlope <-pHSlope%>%
  mutate(pH_insitu = seacarb::pHinsi(pH = pH, ALK = 2200, Tinsi = TEMPINSITU, Tlab = TEMPINLAB, 
                            S = SALINITY,Pt = 0.1, k1k2 = "m10", kf = "dg")) %>%
  select(!pH)%>%
  rename(pH = pH_insitu) %>% # rename it 
  ungroup() %>%
  select(-c(TEMPINLAB, mV, TrisCalDate, TTris, `(Intercept)`, mVTris))

write_csv(pHSlope, here("Data", "Chemistry", "24_hr_pH_data.csv"))

## calculate inflow data using pH slope and flow by each inflow table ##
InflowData <- pHSlope %>%
  filter(TANK_NUM %in% c("Inflow1","Inflow2")) %>%
  select(-c(FLOW_LEFT, FLOW_RIGHT, DO_MG_L, SALINITY, TEMPINSITU))  %>% ### remove the values that I don't need -- You will eventually need to keep TA which is why I dropped these instead of coding for the ones that I need
  rename(pH_inflow = pH,
         TA_inflow = TA) %>%# rename the pH to show that it is inflow pH
  mutate(InflowTable = ifelse(TANK_NUM == "Inflow1",1,2)) %>% # give them inflow numbers to pair easily with the TankID 
  ungroup()%>%
  select(DATE,TIME,InflowTable, pH_inflow, TA_inflow) # drop the Tank ID column to be able to join the data correctly by inflow #

SurfaceArea <- 22.5*22.5 # put in the surface area in cm2 for the bottom of the tank here

## join inflow data, table ID. calculate residence time, delta TA, and NEC ##
Data<-pHSlope %>%
  ungroup()%>%
  filter(!TANK_NUM %in% c("Inflow1","Inflow2"))%>% # filter out the inflow data now
  mutate(TANK_NUM = as.numeric(TANK_NUM))%>% # convert to numeric since the inflow data is now dropped
  left_join(TableID) %>%
  left_join(InflowData) %>% # join with the inflow data for easier calculations of rates
  mutate(DATETIME = ymd_hms(paste(DATE,TIME)), # make a datetime
         pHDiff = pH - pH_inflow, # calculate the difference between the inflow and the pH in each tank 
         totalflow = FLOW_LEFT+FLOW_RIGHT,
         residence_time = (1/totalflow)*(10000/60),# convert ml/min to hours by multiplying by the volumne of water in ml and divide by 60
         deltaTA = TA_inflow - TA, # calculate the difference between in and outflow
         NEC = (deltaTA/2)*(1.025)*(10)*(1/residence_time)*(1/SurfaceArea) ### for a real rate should probably normalize the delta TA to the delta control just like in respo
  )

## pH diffs for EACH tank in 24 hr experiment ##
Data$TREATMENT <- factor(Data$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"), 
                                                    labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"))

tank_pH_diffs_data <- Data %>%
  group_by(DATETIME, TANK_NUM, TREATMENT) %>%
  select(DATETIME, TANK_NUM, TREATMENT, pHDiff) %>%
  summarize(avg_pH_diff = mean(pHDiff, na.rm = TRUE))
  

tank_pH_diffs <- tank_pH_diffs_data %>%
  ggplot(aes(x = DATETIME, y = avg_pH_diff)) +
  geom_rect(aes(xmin = ymd_hms("2024-06-02 12:00:00"), xmax = ymd_hms("2024-06-02 15:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-02 15:00:00"), xmax = ymd_hms("2024-06-02 18:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-02 18:00:00"), xmax = ymd_hms("2024-06-02 21:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-02 21:00:00"), xmax = ymd_hms("2024-06-03 00:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+ 
  geom_rect(aes(xmin = ymd_hms("2024-06-03 00:00:00"), xmax = ymd_hms("2024-06-03 03:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-03 03:00:00"), xmax = ymd_hms("2024-06-03 06:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightgrey", color = NA)+
  geom_rect(aes(xmin = ymd_hms("2024-06-03 06:00:00"), xmax = ymd_hms("2024-06-03 09:00:00"), ymin = -Inf, ymax = Inf),
            alpha = 1/5,
            fill = "lightyellow", color = NA) +
  geom_point(aes(color = TREATMENT)) +
  geom_line(aes(group = TANK_NUM, alpha = 0.5, color = TREATMENT)) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic()+
  annotate("text", x = ymd_hms("2024-06-03 08:00:00"), y = 0.15, label = "Overcast \n Day", size = 4) +
  labs(x= "Date & Time",
       y = "Change in pH") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none") +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"))
tank_pH_diffs
#ggsave(plot = tank_pH_diffs, filename = here("Output", "Supp_Fig_1.png"))

