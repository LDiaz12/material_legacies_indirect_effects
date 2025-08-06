library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(car)
library(lubridate)
library(ggridges)
library(agricolae)
library(lme4)
library(lmerTest)
library(moments)
library(performance)
library(ggpubr)
library(emmeans)
library(agricolae)
library(patchwork)
library(lmtest)
library(nlme)
#install.packages("multcomp")
#library(multcomp)

## bring in pH calibration files and raw data files
pHcalib<-read_csv(here("Data","Chemistry", "TrisCalSummer2024.csv"))
pHData<-read_csv(here("Data", "Chemistry", "CarbonateChemistry.csv"))
TableID<-read_csv(here("Data", "TableID.csv"))
DOC_full_data <- read_csv(here("Data", "DOC", "DOC_full_data.csv"))

pHSlope<-pHcalib %>%
  nest_by(TrisCalDate)%>%
  mutate(fitpH = list(lm(mVTris~TTris, data = pHcalib))) %>% # linear regression of mV and temp of the tris
  reframe(broom::tidy(fitpH)) %>% 
  select(TrisCalDate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>% # put slope and intercept in their own column
  right_join(.,pHData) %>% # join with the pH sample data
  mutate(mVTris = TEMPINLAB*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  drop_na(TEMPINSITU)%>%
  drop_na(mV) %>%
  mutate(pH = pH(Ex=mV,Etris=mVTris,S=SALINITY,T=TEMPINLAB)) %>%  # calculate pH of the samples using the pH seacarb function
  drop_na(TEMPINSITU, TEMPINLAB, SALINITY, pH) %>%
  mutate(pH_insitu = pHinsi(pH = pH, ALK = 2200, Tinsi = TEMPINSITU, Tlab = TEMPINLAB, 
                            S = SALINITY, Pt = 0.1, k1k2 = "m10", kf = "dg")) %>%
  select(!pH)%>%
  rename(pH = pH_insitu) %>% 
  ungroup() %>%
  select(-c(mV, TrisCalDate, TTris, `(Intercept)`, mVTris)) ## warnings are fine, ignore them

#Create new data frame for DIC calculation where you remove NAs from the variables you need 
DIC_Data_calc <- pHSlope %>%
  drop_na(pH, TA)%>%
  mutate(TA_mol_kg = (TA/1e6)) #Convert TA from umol/kg to mol/kg which is what seacarb needs 


#Create the seacarb table. Flag = 8 means you are going to use pH and TA to calculate the rest of the variables
#Seacarb is not tidy so need to use DIC_Data_calc$parameter in order to call them in 
carb_chem_table <- carb(flag = 8, var1 = DIC_Data_calc$pH, var2 = DIC_Data_calc$TA_mol_kg, S = DIC_Data_calc$SALINITY, T = DIC_Data_calc$TEMPINSITU, P = 0, Patm = 0, Pt = 0, Sit = 0, 
                        pHscale = "T", kf = "dg", k1k2 = "m10", ks = "d")

#write_csv(carb_chem_table, here("Data", "Chemistry", "carb_chem_table.csv"))

#Remove just DIC column from carb chem table 
DIC_column <- carb_chem_table %>%
  select(DIC) %>%
  mutate(DIC_µmol_kg = (DIC*1e6)) 

#Join the DIC column back into the DIC_Data_calc data frame 
DIC_Data_calc2 <- DIC_Data_calc %>%
  bind_cols(DIC_column)

#Rejoin the newly calculated DIC (mol/kg) data frame back to the main 'Data' data frame
pHSlope2 <- pHSlope %>%
  left_join(DIC_Data_calc2)%>% 
  select(-c(TA_mol_kg,DIC)) %>% # now we can calculate NEP 
  mutate(TREATMENT = ifelse(is.na(TREATMENT),"Inflow",TREATMENT)) %>% # rename inflow treatment
  full_join(DOC_full_data %>%
              select(!DATETIME)) %>%
  select(-c(NPOC_mg_L, TN_mg_L))


# remove the inflow data and join it with the tanks that had that specific inflow water
InflowData <- pHSlope2 %>%
  filter(TANK_NUM %in% c("Inflow1","Inflow2")) %>%
  select(-c(FLOW_LEFT, FLOW_RIGHT, Notes, DO_MG_L, SALINITY, TEMPINSITU))  %>%
  rename(pH_inflow = pH,
         TA_inflow = TA,
         DIC_inflow = DIC_µmol_kg,
         DOC_inflow = NPOC_uM) %>%
  mutate(INFLOW_TABLE = ifelse(TANK_NUM == "Inflow1",1,2)) %>% # give them inflow numbers to pair easily with the TankID 
  ungroup()%>%
  select(DATE,TIME, INFLOW_TABLE, pH_inflow, TA_inflow, DIC_inflow, DOC_inflow) # drop the Tank ID column to be able to join the data correctly by inflow #

#InflowData <- InflowData %>%
 # mutate(pH_inflow = ifelse(DATE == "20240606" & INFLOW_TABLE == "2" & pH_inflow > 8.05, NA, pH_inflow),
  #       TA_inflow = ifelse(DATE == "20240606" & INFLOW_TABLE == "1" & TA_inflow < 2300, NA, TA_inflow),
   #      DOC_inflow = ifelse(DATE == "20240606" & INFLOW_TABLE == "1", NA, DOC_inflow), 
    #     DOC_inflow = ifelse(DATE == "20240611" & INFLOW_TABLE == "1" & DOC_inflow > 300, NA, DOC_inflow))


SurfaceArea <- 22.5*22.5 #cm^2

Data<-pHSlope2 %>%
  ungroup()%>%
  filter(!TANK_NUM %in% c("Inflow1","Inflow2"))%>% # filter out the inflow data
  mutate(TANK_NUM = as.numeric(TANK_NUM))%>% # convert to numeric since the inflow data is now dropped
  left_join(TableID) %>%
  left_join(InflowData) %>% # join with the inflow data for easier calculations of rates
  mutate(DATETIME = ymd_hms(paste(DATE,TIME)), # make a datetime
         deltapH = pH - pH_inflow, # calculate the difference between the inflow and the pH in each tank 
         deltaDOC = NPOC_uM - DOC_inflow, # calculate difference from DOC inflow and each tank
         totalflow = FLOW_RIGHT+FLOW_LEFT,
         residence_time = (1/totalflow)*(10000/60),# convert ml/min to hours by multiplying by the volume of water in ml (10L tank; 10,000mL) and divide by 60 mins
         flowrate = (totalflow/60), # convert mL/min to mL/sec
         deltaTA = TA_inflow - TA, # calculate the difference between in and outflow
         deltaDIC = DIC_inflow - DIC_µmol_kg, # calculate delta DIC to calculate NEP 
         NEC = (deltaTA/2)*(1.025)*(10)*(1/residence_time)*(1/SurfaceArea), ### for a real rate should probably normalize the delta TA to the delta control just like in respo
         NEP = ((deltaDIC)*(1.025)*(10)*(1/residence_time)*(1/SurfaceArea)) - NEC)


# table of light and temp per day # 
Data1 <- Data %>%
  filter(!TEMPINSITU > 36)

light_per_day <- Data1 %>%
  group_by(DATE) %>%
  summarize(avg_light = mean(LIGHT_NM, na.rm = TRUE))
light_per_day

temp_per_day <- Data1 %>%
  group_by(DATE) %>%
  summarize(avg_temp = mean(TEMPINSITU, na.rm = TRUE))
temp_per_day


temp_per_day_plot <- light_temp_per_day %>%
  ggplot(aes(x = DATETIME, y = avg_temp)) + 
  geom_point()
temp_per_day_plot

temp_per_day_plot <- light_temp_per_day %>%
  ggplot(aes(x = DATETIME, y = avg_light)) + 
  geom_point()
temp_per_day_plot

## write full chemistry data file ##
# temp, flow, salinity, pH, light, DIC, pH inflow, TA inflow, DIC inflow, 
# delta pH, delta TA, delta DIC, NEC, NEP, DOC 


#write_csv(Data, here("Data", "Chemistry", "Full_Carb_Chem_Data.csv"))
full_carb_chem_data <- read_csv(here("Data", "Chemistry", "Full_Carb_Chem_Data.csv"))
full_carb_chem_data$INFLOW_TABLE <- as.factor(full_carb_chem_data$INFLOW_TABLE)

#checking for outliers in inflow data per inflow table for TA, pH, and DOC 
ggplot(full_carb_chem_data) + 
  geom_point(aes(x = DATETIME, y = DOC_inflow, color = INFLOW_TABLE)) + 
  facet_wrap(~DATE)
# 06/06/2024 looks like it was a bad date for each of the parameters. it was raining all day and overcast

calc_avg_res_flow <- Data %>%
  group_by(TANK_NUM) %>%
  summarize(avg_res_time = mean(residence_time, na.rm = TRUE),
            avg_flow_rate = mean(flowrate, na.rm = TRUE)) %>%
  drop_na()
calc_avg_res_flow

avg_total_flow <- Data %>%
  summarize(avg_total_flow = mean(flowrate, na.rm = TRUE),
            flow_error = sd(flowrate, na.rm = TRUE)/sqrt(n()))
avg_total_flow


## create one summary data frame with means and ranges for TA, pH, and DOC ##
# start by isolating 12:00 and 21:00 time periods #
## filtering for 12:00 and 21:00 sampling times and reframing to add daily mean and daily range between 12 and 9 

chem_reframe <- Data %>% 
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  select(DATE, TIME, TANK_NUM, TREATMENT, TA, pH, DIC_µmol_kg, NPOC_uM, deltapH, deltaTA, deltaDOC, 
         deltaDIC, NEC, NEP) %>%
  reframe(TA_range = TA[TIME == hms("12:00:00")] - TA[TIME == hms("21:00:00")],
          TA_dailymean = mean(TA, na.rm = TRUE),
          deltaTA_range = deltaTA[TIME == hms("12:00:00")] - deltaTA[TIME == hms("21:00:00")],
          deltaTA_dailymean = mean(deltaTA, na.rm = TRUE),
          pH_range = pH[TIME == hms("12:00:00")] - pH[TIME == hms("21:00:00")],
          pH_dailymean = mean(pH, na.rm = TRUE),
          deltapH_range = deltapH[TIME == hms("12:00:00")] - deltapH[TIME == hms("21:00:00")],
          deltapH_dailymean = mean(deltapH, na.rm = TRUE),
          DOC_range = NPOC_uM[TIME == hms("12:00:00")] - NPOC_uM[TIME == hms("21:00:00")],
          DOC_dailymean = mean(NPOC_uM, na.rm = TRUE),
          deltaDOC_range = deltaDOC[TIME == hms("12:00:00")] - deltaDOC[TIME == hms("21:00:00")],
          deltaDOC_dailymean = mean(deltaDOC, na.rm = TRUE),
          NEC_range = NEC[TIME == hms("12:00:00")] - NEC[TIME == hms("21:00:00")],
          NEC_dailymean = mean(NEC, na.rm = TRUE),
          NEP_range = NEP[TIME == hms("12:00:00")] - NEP[TIME == hms("21:00:00")],
          NEP_dailymean = mean(NEP, na.rm = TRUE))

## reframe chem data only during DAY time 
chem_reframe_DAY <- Data %>%
  filter(TIME %in% "12:00:00") %>%
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  select(DATE, TIME, TANK_NUM, TREATMENT, TA, pH, DIC_µmol_kg, NPOC_uM, deltapH, deltaTA, deltaDOC, 
         deltaDIC, NEC, NEP) %>%
  rename(DOC = NPOC_uM)
## review chem reframe DAY and look for outliers 
chem_reframe_DAY_plots <- chem_reframe_DAY %>%
  select(DATE, TIME, TREATMENT, TANK_NUM, TA:NEP) %>%
  pivot_longer(cols = TA:NEP) %>%
  ggplot(aes(x = TREATMENT, y = value)) +
  geom_point() + 
  facet_wrap(~name, scales = "free")
chem_reframe_DAY_plots
# clean up outliers 
chem_reframe_DAY_clean <- chem_reframe_DAY %>%
  mutate(deltaDIC = ifelse(deltaDIC < -250, NA, deltaDIC),
         deltaDOC = ifelse(deltaDOC > 200, NA, deltaDOC),
         deltaTA = ifelse(deltaTA < -400, NA, deltaTA), 
         DIC_µmol_kg = ifelse(DIC_µmol_kg > 2250, NA, DIC_µmol_kg), 
         DOC = ifelse(DOC > 350, NA, DOC), 
         NEC = ifelse(NEC < -2, NA, NEC), 
         NEC = ifelse(TREATMENT == "Rubble_Dom" & NEC > 2, NA, NEC), 
         NEP = ifelse(NEP < -2.1, NA, NEP),
         NEP = ifelse(NEP > 5, NA, NEP), 
         pH = ifelse(pH < 7.8, NA, pH), 
         TA = ifelse(TA > 2700, NA, TA))

# reframe chem data only during NIGHT 
chem_reframe_NIGHT <- Data %>%
  filter(TIME %in% "21:00:00") %>%
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  select(DATE, TIME, TANK_NUM, TREATMENT, TA, pH, DIC_µmol_kg, NPOC_uM, deltapH, deltaTA, deltaDOC, 
         deltaDIC, NEC, NEP) %>%
  rename(DOC = NPOC_uM)
# review chem reframe NIGHT for outliers and clean 
chem_reframe_NIGHT_plots <- chem_reframe_NIGHT %>%
  select(DATE, TIME, TREATMENT, TANK_NUM, TA:NEP) %>%
  pivot_longer(cols = TA:NEP) %>%
  ggplot(aes(x = TREATMENT, y = value)) +
  geom_point() + 
  facet_wrap(~name, scales = "free")
chem_reframe_NIGHT_plots
# clean up outliers 
chem_reframe_NIGHT_clean <- chem_reframe_NIGHT %>%
  mutate(deltaDOC = ifelse(deltaDOC > 300, NA, deltaDOC),
         deltaTA = ifelse(deltaTA > 100, NA, deltaTA), 
         DIC_µmol_kg = ifelse(DIC_µmol_kg < 1900, NA, DIC_µmol_kg), 
         DOC = ifelse(DOC > 400, NA, DOC), 
         NEC = ifelse(NEC > 1, NA, NEC), 
         NEP = ifelse(NEP > 2, NA, NEP), 
         TA = ifelse(TA < 2250, NA, TA))

# make summary data of tanks at night # 
chem_summary_data_NIGHT <- chem_reframe_NIGHT_clean %>%
  group_by(TREATMENT, TANK_NUM) %>%
  summarize(TA_mean = mean(TA, na.rm = TRUE),
            TA_se = sd(TA, na.rm = TRUE)/sqrt(n()),
            deltaTA_mean = mean(deltaTA, na.rm = TRUE),
            deltaTA_se = sd(deltaTA, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH, na.rm = TRUE),
            pH_se = sd(pH, na.rm = TRUE)/sqrt(n()),
            deltapH_mean = mean(deltapH, na.rm = TRUE),
            deltapH_se = sd(deltapH, na.rm = TRUE)/sqrt(n()),
            DOC_mean = mean(DOC, na.rm = TRUE),
            DOC_se = sd(DOC, na.rm = TRUE)/sqrt(n()),
            deltaDOC_mean = mean(deltaDOC, na.rm = TRUE),
            deltaDOC_se = sd(deltaDOC, na.rm = TRUE)/sqrt(n()),
            NEC_mean = mean(NEC, na.rm = TRUE),
            NEC_se = sd(NEC, na.rm = TRUE)/sqrt(n()),
            NEP_mean = mean(NEP, na.rm = TRUE),
            NEP_se = sd(NEP, na.rm = TRUE)/sqrt(n()))
#write_csv(chem_summary_data_NIGHT, here("Data", "Chemistry", "chem_summary_data_NIGHT.csv"))

# make summary data of tanks during the day # 
chem_summary_data_DAY <- chem_reframe_DAY_clean %>%
  group_by(TREATMENT, TANK_NUM) %>%
  summarize(TA_mean = mean(TA, na.rm = TRUE),
            TA_se = sd(TA, na.rm = TRUE)/sqrt(n()),
            deltaTA_mean = mean(deltaTA, na.rm = TRUE),
            deltaTA_se = sd(deltaTA, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH, na.rm = TRUE),
            pH_se = sd(pH, na.rm = TRUE)/sqrt(n()),
            deltapH_mean = mean(deltapH, na.rm = TRUE),
            deltapH_se = sd(deltapH, na.rm = TRUE)/sqrt(n()),
            DOC_mean = mean(DOC, na.rm = TRUE),
            DOC_se = sd(DOC, na.rm = TRUE)/sqrt(n()),
            deltaDOC_mean = mean(deltaDOC, na.rm = TRUE),
            deltaDOC_se = sd(deltaDOC, na.rm = TRUE)/sqrt(n()),
            NEC_mean = mean(NEC, na.rm = TRUE),
            NEC_se = sd(NEC, na.rm = TRUE)/sqrt(n()),
            NEP_mean = mean(NEP, na.rm = TRUE),
            NEP_se = sd(NEP, na.rm = TRUE)/sqrt(n()))
#write_csv(chem_summary_data_DAY, here("Data", "Chemistry", "chem_summary_data_DAY.csv"))

chem_reframe_plots <- chem_reframe %>%
  select(TREATMENT, TANK_NUM, TA_range:NEP_dailymean) %>%
  pivot_longer(cols = TA_range:NEP_dailymean) %>%
  ggplot(aes(x = TREATMENT, y = value)) +
  geom_point() + 
  facet_wrap(~name, scales = "free")
chem_reframe_plots

chem_reframe_clean <- chem_reframe %>%
  mutate(deltaTA_dailymean = ifelse(deltaTA_dailymean < -250, NA, deltaTA_dailymean),
         deltapH_dailymean = ifelse(deltapH_dailymean > 0.4, NA, deltapH_dailymean), 
         deltapH_range = ifelse(deltapH_range < -0.15, NA, deltapH_range), 
         deltaDOC_dailymean = ifelse(deltaDOC_dailymean > 300, NA, deltaDOC_dailymean),
         deltaDOC_range = ifelse(deltaDOC_range > 250, NA, deltaDOC_range),
         DOC_dailymean = ifelse(DOC_dailymean > 400, NA, DOC_dailymean),
         DOC_range = ifelse(DOC_range < -200, NA, DOC_range),
         NEC_dailymean = ifelse(TREATMENT ==  "Coral_Dom"&NEC_dailymean < -9, NA, NEC_dailymean),
         NEC_range = ifelse(NEC_range < -1, NA, NEC_range),
         NEP_range = ifelse(NEP_range < -2, NA, NEP_range), 
         NEP_dailymean = ifelse(NEP_dailymean < -5, NA, NEP_dailymean),
         TA_dailymean = ifelse(TA_dailymean > 2500, NA, TA_dailymean), 
         TA_dailymean = ifelse(TREATMENT == "Coral_Dom"&TA_dailymean < 2200, NA, TA_dailymean),
         TA_range = ifelse(TA_range > 400, NA, TA_range))

chem_reframe_clean %>%
  select(TREATMENT, TANK_NUM, TA_range:NEP_dailymean) %>%
  pivot_longer(cols = TA_range:NEP_dailymean) %>%
  ggplot(aes(x = TREATMENT, y = value)) +
  geom_point() + 
  facet_wrap(~name, scales = "free")

chem_reframe_clean_ALGAE <- chem_reframe_clean %>%
  filter(TREATMENT == "Algae_Dom") %>%
  select(TANK_NUM, TA_range:NEP_dailymean)
chem_reframe_clean_ALGAE

chem_reframe_clean_CORAL <- chem_reframe_clean %>%
  filter(TREATMENT == "Coral_Dom") %>%
  select(TANK_NUM, TA_range:NEP_dailymean)
chem_reframe_clean_CORAL

chem_reframe_clean_RUBBLECCA <- chem_reframe_clean %>%
  filter(TREATMENT == "Rubble_Dom") %>%
  select(TANK_NUM, TA_range:NEP_dailymean)
chem_reframe_clean_RUBBLECCA

chem_reframe_clean_CONTROL <- chem_reframe_clean %>%
  filter(TREATMENT == "Control") %>%
  select(TANK_NUM, TA_range:NEP_dailymean)
chem_reframe_clean_CONTROL

chem_summary_data <- chem_reframe_clean %>%
  group_by(TREATMENT, TANK_NUM) %>%
  summarize(TA_rangemean = mean(TA_range, na.rm = TRUE), 
            TA_rangese = sd(TA_range, na.rm = TRUE)/sqrt(n()),
            TA_mean = mean(TA_dailymean, na.rm = TRUE),
            TA_se = sd(TA_dailymean, na.rm = TRUE)/sqrt(n()),
            deltaTA_rangemean = mean(deltaTA_range, na.rm = TRUE), 
            deltaTA_rangese = sd(deltaTA_range, na.rm = TRUE)/sqrt(n()),
            deltaTA_mean = mean(deltaTA_dailymean, na.rm = TRUE),
            deltaTA_se = sd(deltaTA_dailymean, na.rm = TRUE)/sqrt(n()),
            pH_rangemean = mean(pH_range, na.rm = TRUE), 
            pH_rangese = sd(pH_range, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH_dailymean, na.rm = TRUE),
            pH_se = sd(pH_dailymean, na.rm = TRUE)/sqrt(n()),
            deltapH_rangemean = mean(deltapH_range, na.rm = TRUE), 
            deltapH_rangese = sd(deltapH_range, na.rm = TRUE)/sqrt(n()),
            deltapH_mean = mean(deltapH_dailymean, na.rm = TRUE),
            deltapH_se = sd(deltapH_dailymean, na.rm = TRUE)/sqrt(n()),
            DOC_rangemean = mean(DOC_range, na.rm = TRUE), 
            DOC_rangese = sd(DOC_range, na.rm = TRUE)/sqrt(n()),
            DOC_mean = mean(DOC_dailymean, na.rm = TRUE),
            DOC_se = sd(DOC_dailymean, na.rm = TRUE)/sqrt(n()),
            deltaDOC_rangemean = mean(deltaDOC_range, na.rm = TRUE), 
            deltaDOC_rangese = sd(deltaDOC_range, na.rm = TRUE)/sqrt(n()),
            deltaDOC_mean = mean(deltaDOC_dailymean, na.rm = TRUE),
            deltaDOC_se = sd(deltaDOC_dailymean, na.rm = TRUE)/sqrt(n()),
            NEC_rangemean = mean(NEC_range, na.rm = TRUE),
            NEC_rangese = sd(NEC_range, na.rm = TRUE)/sqrt(n()),
            NEC_mean = mean(NEC_dailymean, na.rm = TRUE),
            NEC_se = sd(NEC_dailymean, na.rm = TRUE)/sqrt(n()),
            NEP_rangemean = mean(NEP_range, na.rm = TRUE), 
            NEP_rangese = sd(NEP_range, na.rm = TRUE)/sqrt(n()),
            NEP_mean = mean(NEP_dailymean, na.rm = TRUE),
            NEP_se = sd(NEP_dailymean, na.rm = TRUE)/sqrt(n()))

chem_summary_data_2 <- chem_reframe_date_removed_clean %>%
  group_by(TREATMENT, TANK_NUM) %>%
  summarize(TA_rangemean = mean(TA_range, na.rm = TRUE), 
            TA_rangese = sd(TA_range, na.rm = TRUE)/sqrt(n()),
            TA_mean = mean(TA_dailymean, na.rm = TRUE),
            TA_se = sd(TA_dailymean, na.rm = TRUE)/sqrt(n()),
            deltaTA_rangemean = mean(deltaTA_range, na.rm = TRUE), 
            deltaTA_rangese = sd(deltaTA_range, na.rm = TRUE)/sqrt(n()),
            deltaTA_mean = mean(deltaTA_dailymean, na.rm = TRUE),
            deltaTA_se = sd(deltaTA_dailymean, na.rm = TRUE)/sqrt(n()),
            pH_rangemean = mean(pH_range, na.rm = TRUE), 
            pH_rangese = sd(pH_range, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH_dailymean, na.rm = TRUE),
            pH_se = sd(pH_dailymean, na.rm = TRUE)/sqrt(n()),
            deltapH_rangemean = mean(deltapH_range, na.rm = TRUE), 
            deltapH_rangese = sd(deltapH_range, na.rm = TRUE)/sqrt(n()),
            deltapH_mean = mean(deltapH_dailymean, na.rm = TRUE),
            deltapH_se = sd(deltapH_dailymean, na.rm = TRUE)/sqrt(n()),
            DOC_rangemean = mean(DOC_range, na.rm = TRUE), 
            DOC_rangese = sd(DOC_range, na.rm = TRUE)/sqrt(n()),
            DOC_mean = mean(DOC_dailymean, na.rm = TRUE),
            DOC_se = sd(DOC_dailymean, na.rm = TRUE)/sqrt(n()),
            deltaDOC_rangemean = mean(deltaDOC_range, na.rm = TRUE), 
            deltaDOC_rangese = sd(deltaDOC_range, na.rm = TRUE)/sqrt(n()),
            deltaDOC_mean = mean(deltaDOC_dailymean, na.rm = TRUE),
            deltaDOC_se = sd(deltaDOC_dailymean, na.rm = TRUE)/sqrt(n()),
            NEC_rangemean = mean(NEC_range, na.rm = TRUE),
            NEC_rangese = sd(NEC_range, na.rm = TRUE)/sqrt(n()),
            NEC_mean = mean(NEC_dailymean, na.rm = TRUE),
            NEC_se = sd(NEC_dailymean, na.rm = TRUE)/sqrt(n()),
            NEP_rangemean = mean(NEP_range, na.rm = TRUE), 
            NEP_rangese = sd(NEP_range, na.rm = TRUE)/sqrt(n()),
            NEP_mean = mean(NEP_dailymean, na.rm = TRUE),
            NEP_se = sd(NEP_dailymean, na.rm = TRUE)/sqrt(n()))


            
#write_csv(chem_summary_data, here("Data", "Chemistry", "chem_summary_data.csv"))
#write_csv(chem_summary_data_2, here("Data", "Chemistry", "chem_summary_data_dateremoved.csv"))

tank_chem_means <- chem_reframe_clean %>%
  group_by(TANK_NUM) %>%
  summarize(avg_TA = mean(TA_dailymean, na.rm = TRUE),
            TA_error = sd(TA_dailymean, na.rm = TRUE)/sqrt(n()),
            avg_range_TA = mean(TA_range, na.rm = TRUE),
            TA_range_error = sd(TA_range, na.rm = TRUE)/sqrt(n()),
            avg_pH = mean(pH_dailymean, na.rm = TRUE), 
            pH_error = sd(pH_dailymean, na.rm = TRUE)/sqrt(n()),
            avg_range_pH = mean(pH_range, na.rm = TRUE), 
            pH_range_error = sd(pH_range, na.rm = TRUE)/sqrt(n()),
            avg_DOC = mean(DOC_dailymean, na.rm = TRUE), 
            DOC_error = sd(DOC_dailymean, na.rm = TRUE)/sqrt(n()))
tank_chem_means

######################
## create rect intervals for the NEC plot below for light vs dark times ##
rect_intervals <- tibble::tibble(
  xmin = as.POSIXct(c("2024-06-02 06:00:00", "2024-06-02 18:00:00", 
                      "2024-06-03 06:00:00", "2024-06-03 18:00:00", 
                      "2024-06-04 06:00:00", "2024-06-04 18:00:00",
                      "2024-06-05 06:00:00", "2024-06-05 18:00:00",
                      "2024-06-06 06:00:00", "2024-06-06 18:00:00",
                      "2024-06-07 06:00:00", "2024-06-07 18:00:00",
                      "2024-06-08 06:00:00", "2024-06-08 18:00:00",
                      "2024-06-09 06:00:00", "2024-06-09 18:00:00",
                      "2024-06-10 06:00:00", "2024-06-10 18:00:00",
                      "2024-06-11 06:00:00", "2024-06-11 18:00:00",
                      "2024-06-12 06:00:00", "2024-06-12 18:00:00",
                      "2024-06-13 06:00:00", "2024-06-13 18:00:00",
                      "2024-06-14 06:00:00", "2024-06-14 18:00:00",
                      "2024-06-15 06:00:00", "2024-06-15 18:00:00",
                      "2024-06-16 06:00:00", "2024-06-16 18:00:00",
                      "2024-06-17 06:00:00", "2024-06-17 18:00:00",
                      "2024-06-18 06:00:00", "2024-06-18 18:00:00",
                      "2024-06-19 06:00:00", "2024-06-19 18:00:00",
                      "2024-06-20 06:00:00", "2024-06-20 18:00:00",
                      "2024-06-21 06:00:00", "2024-06-21 18:00:00",
                      "2024-06-22 06:00:00", "2024-06-22 18:00:00",
                      "2024-06-23 06:00:00", "2024-06-23 18:00:00",
                      "2024-06-24 06:00:00", "2024-06-24 18:00:00",
                      "2024-06-25 06:00:00", "2024-06-25 18:00:00",
                      "2024-06-26 06:00:00", "2024-06-26 18:00:00")),
  xmax = as.POSIXct(c("2024-06-02 18:00:00", "2024-06-03 06:00:00", 
                      "2024-06-03 18:00:00", "2024-06-04 06:00:00", 
                      "2024-06-04 18:00:00", "2024-06-05 06:00:00",
                      "2024-06-05 18:00:00", "2024-06-06 06:00:00",
                      "2024-06-06 18:00:00", "2024-06-07 06:00:00",
                      "2024-06-07 18:00:00", "2024-06-08 06:00:00",
                      "2024-06-08 18:00:00", "2024-06-09 06:00:00",
                      "2024-06-09 18:00:00", "2024-06-10 06:00:00",
                      "2024-06-10 18:00:00", "2024-06-11 06:00:00",
                      "2024-06-11 18:00:00", "2024-06-12 06:00:00",
                      "2024-06-12 18:00:00", "2024-06-13 06:00:00",
                      "2024-06-13 18:00:00", "2024-06-14 06:00:00",
                      "2024-06-14 18:00:00", "2024-06-15 06:00:00",
                      "2024-06-15 18:00:00", "2024-06-16 06:00:00",
                      "2024-06-16 18:00:00", "2024-06-17 06:00:00",
                      "2024-06-17 18:00:00", "2024-06-18 06:00:00",
                      "2024-06-18 18:00:00", "2024-06-19 06:00:00",
                      "2024-06-19 18:00:00", "2024-06-20 06:00:00",
                      "2024-06-20 18:00:00", "2024-06-21 06:00:00",
                      "2024-06-21 18:00:00", "2024-06-22 06:00:00",
                      "2024-06-22 18:00:00", "2024-06-23 06:00:00",
                      "2024-06-23 18:00:00", "2024-06-24 06:00:00",
                      "2024-06-24 18:00:00", "2024-06-25 06:00:00",
                      "2024-06-25 18:00:00", "2024-06-26 06:00:00",
                      "2024-06-26 18:00:00", "2024-06-27 06:00:00")),
  fill = c("lightyellow", "lightgrey", "lightyellow", "lightgrey", "lightyellow",
           "lightgrey", "lightyellow", "lightgrey", "lightyellow", "lightgrey", 
           "lightyellow", "lightgrey", "lightyellow","lightgrey", "lightyellow", 
           "lightgrey", "lightyellow","lightgrey", "lightyellow", "lightgrey", 
           "lightyellow","lightgrey", "lightyellow", "lightgrey", "lightyellow",
           "lightgrey", "lightyellow", "lightgrey", "lightyellow","lightgrey", 
           "lightyellow", "lightgrey", "lightyellow","lightgrey", "lightyellow", 
           "lightgrey", "lightyellow","lightgrey", "lightyellow", "lightgrey", 
           "lightyellow", "lightgrey", "lightyellow", "lightgrey", "lightyellow",
           "lightgrey", "lightyellow","lightgrey", "lightyellow", "lightgrey"))

### Plot raw TA data ###
raw_TA_plot <- Data %>%
  group_by(TREATMENT, DATETIME)%>%
  ggplot(aes(x = DATETIME, y = TA, color = TREATMENT, na.rm = TRUE)) +
  geom_rect(data = rect_intervals, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 1/5, color = NA, inherit.aes = FALSE) +
  scale_x_datetime(date_labels = "%D\n%T", breaks= seq(min(Data$DATETIME), max(Data$DATETIME), 
                                                             by = "3 days")) + # is there a way to change the time from 5:00 to 12:00?
  scale_fill_identity() +
  geom_point() +
  geom_jitter(data = Data, aes(x = DATETIME, y = TA), alpha = 0.7) +
  theme_classic() +
  labs(x="Date & Time",
       y = "Total Alkalinity (umol/kg)") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
raw_TA_plot
#ggsave(plot = raw_TA_plot, filename = here("Output", "TA_NECPlots", "raw_TA_plot.png"), width = 9, height = 9)

### Plot mean TA ###
TA_plot <- Data %>%
  group_by(TREATMENT, DATETIME)%>%
  summarise(mean_TA = mean(TA, na.rm = TRUE),
            se_TA = sd(TA, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = DATETIME, y = mean_TA, color = TREATMENT, na.rm = TRUE)) +
  geom_rect(data = rect_intervals, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 1/5, color = NA, inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(x="Date & Time",
       y = "Mean Total Alkalinity (umol/kg)") +
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
TA_plot

## TA RANGE ## 
# reorder treatments in order I want presented on graphs by changing factor levels # 
chem_summary_data$TREATMENT <- factor(chem_summary_data$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

TA_range_plot <- chem_summary_data %>%
  ggplot(aes(x = TREATMENT, y = TA_rangemean, color = TREATMENT)) +
  geom_point() +
  labs(x = "Treatment", y = expression(bold("Daily Mean Total Alkalinity Range" ~ (µmol ~ kg^-1)))) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1, color = "black") +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"))
TA_range_plot
#ggsave(plot = TA_range_plot, filename = here("Output", "TA_NECPlots", "TA_range_plot.png"), width = 9, height = 7)

# range TA stats #
TA_range_model <- lmer(TA_range ~ TREATMENT + (1|TANK_NUM), data = chem_reframe_clean)
check_model(TA_range_model)  
summary(TA_range_model)
anova(TA_range_model) # significant effect of treatment on range in TA 

## TA MEAN ## 
# mean TA plot #
TA_mean_plot <- chem_summary_data %>%
  ggplot(aes(x = TREATMENT, y = TA_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daily Mean TA" ~ (µmol ~ kg^-1)))) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1, color = "black") +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"))
TA_mean_plot
#ggsave(plot = TA_mean_plot, filename = here("Output", "TA_mean_plot.png"), width = 9, height = 6)

## TA daily mean stats ##
TA_mean_model <- lmer(TA_dailymean ~ TREATMENT + (1|TANK_NUM), data=TA_data2_filtered)
check_model(TA_mean_model)
summary(TA_mean_model) # rubble dom and coral dom community significant at p < 0.05 for both
anova(TA_mean_model) # significant effect of benthic community on mean TA throughout experimental period 
# this makes sense given the high calcification rates of corals and CCA

# model without random effects for post hoc groupings
TA_mean_model2 <- lm(TA_dailymean ~ TREATMENT, data=TA_data2_filtered)
HSD.test(TA_mean_model2, "TREATMENT", console = TRUE)
# coral and rubble dom communities are most similar, algae is most similar to controls which makes 
# sense bc turbs are non calcifying 

###################################
### NEC DATA ### 

## DAILY MEAN NEC ## 
chem_reframe_clean$TREATMENT <- factor(chem_reframe_clean$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

## NEC DAY vs NEC NIGHT ##
# NEC day

NEC_data_day <- chem_reframe_DAY_clean %>% # plot NEC data during the DAY using cleaned data (cleaned 6/6/25)
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  select(DATE, TIME, TANK_NUM, TREATMENT, NEC) %>%
  reframe(NEC_day_mean = mean(NEC, na.rm = TRUE)) %>%
  drop_na()

NEC_day_treatment_means <- NEC_data_day %>%
  group_by(TREATMENT) %>%
  summarize(mean_NEC_treatment = mean(NEC_day_mean, na.rm = TRUE),
            se_NEC_treatment = sd(NEC_day_mean, na.rm = TRUE)/sqrt(n()))

NEC_data_day$TREATMENT <- factor(NEC_data_day$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEC_day_mean_plot <- NEC_data_day %>%
  ggplot(aes(x = TREATMENT, y = NEC_day_mean, color = TREATMENT)) +
  geom_hline(yintercept = 0, lty = 2) +
  #ylim(c(-1.5,1.5)) +
  labs(x = "", y = expression(bold("NCC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1)))) +
  ggtitle("Daytime Mean") +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = NEC_data_day, aes(x = TREATMENT, y = NEC_day_mean), alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  #geom_text(data = NEC_day_treatment_means, 
            #aes(x = TREATMENT, y = mean_NEC_treatment, 
                #label = paste0("", round(mean_NEC_treatment, 2), "±", round(se_NEC_treatment, 2))),
            #vjust = -3, hjust = 0.2, color = "black", size = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
NEC_day_mean_plot
#ggsave(plot = NEC_day_mean_plot, filename = here("Output", "TA_NECPlots", "NEC_day_mean_plot.png"), width = 9, height = 9)

NEC_daytime_model <- lmer(NEC_day_mean ~ TREATMENT + (1|TANK_NUM), data=NEC_data_day)
check_model(NEC_daytime_model)
summary(NEC_daytime_model)
anova(NEC_daytime_model) # sig effect of community on daytime NEC rates. p = 0.0001494

emmeans(NEC_daytime_model, pairwise ~ "TREATMENT", adjust = "Tukey")
NEC_daytime_model_noRandom <- lm(NEC_day_mean ~ TREATMENT, data = NEC_data_day)
HSD.test(NEC_daytime_model_noRandom, "TREATMENT", console=TRUE)

# NEC night
NEC_data_night <- chem_reframe_NIGHT_clean %>% # plot NEC data during the NIGHT using cleaned data (cleaned 6/6/25)
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  select(DATE, TIME, TANK_NUM, TREATMENT, NEC) %>%
  reframe(NEC_night_mean = mean(NEC, na.rm = TRUE)) %>%
  drop_na()

NEC_night_treatment_means <- NEC_data_night %>%
  group_by(TREATMENT) %>%
  summarize(mean_NEC_treatment = mean(NEC_night_mean, na.rm = TRUE),
            se_NEC_treatment = sd(NEC_night_mean, na.rm = TRUE)/sqrt(n()))

NEC_data_night$TREATMENT <- factor(NEC_data_night$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEC_night_mean_plot <- NEC_data_night %>%
  ggplot(aes(x = TREATMENT, y = NEC_night_mean, color = TREATMENT)) +
  geom_hline(yintercept = 0, lty=2) +
  ggtitle("Nighttime Mean") +
  #ylim(c(-0.5,0.5)) +
  labs(x = "", y = "") +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = NEC_data_night, aes(x = TREATMENT, y = NEC_night_mean), alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  #geom_text(data = NEC_night_treatment_means, 
  #          aes(x = TREATMENT, y = mean_NEC_treatment, 
   #             label = paste0("", round(mean_NEC_treatment, 2), "±", round(se_NEC_treatment, 2))),
    #        vjust = -3, hjust = 0.2, color = "black", size = 4) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none")
NEC_night_mean_plot
#ggsave(plot = NEC_night_mean_plot, filename = here("Output", "TA_NECPlots", "NEC_night_mean_plot.png"), width = 9, height = 9)

NEC_data_night$TREATMENT <- factor(NEC_data_night$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEC_nighttime_model <- lmer(NEC_night_mean ~ TREATMENT + (1|TANK_NUM), data=NEC_data_night)
check_model(NEC_nighttime_model)
summary(NEC_nighttime_model)
anova(NEC_nighttime_model) # sig effect of treatments on nighttime nec. p = 0.007543

emmeans(NEC_nighttime_model, pairwise ~ "TREATMENT", adjust = "Tukey")
NEC_nighttime_model_noRandom <- lm(NEC_night_mean ~ TREATMENT, data = NEC_data_night)
HSD.test(NEC_nighttime_model_noRandom, "TREATMENT", console=TRUE)

##### NEP DATA ANALYSIS ##### 
NEP_data <- Data %>%
  select(DATETIME, DATE, TIME, NEP, TREATMENT, TANK_NUM, INFLOW_TABLE) %>%
  group_by(DATETIME, DATE, TIME, TREATMENT, TANK_NUM, INFLOW_TABLE) %>%
  drop_na()

NEP_data2 <- NEP_data %>% 
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(NEP_range = NEP[TIME == hms("12:00:00")] - NEP[TIME == hms("21:00:00")],
          NEP_dailymean = mean(NEP, na.rm = TRUE)) %>%
  filter(NEP_range > -5) %>%
  filter(!(TREATMENT == "Coral_Dom" & NEP_range < -2)) %>%
  filter(!NEP_range < -2)
#############################
#create plot data

NEP_range_plot <- NEP_plotdata %>%
  ggplot(aes(x = TREATMENT, y = NEP_rangemean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daily Mean NEP Range" ~ (mmol ~ C ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble/CCA-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = NEP_data2, aes(x = TREATMENT, y = NEP_range), alpha = 0.7) +
  geom_errorbar(aes(ymin = NEP_rangemean - NEP_rangese,
                   ymax = NEP_rangemean + NEP_rangese), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEP_range_plot
#ggsave(plot = NEP_range_plot, filename = here("Output", "NEP_Plots", "NEP_range_plot.png"), width = 9, height = 9)

# NEP range stats #
NEP_data2$TREATMENT <- factor(NEP_data2$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))
# reordering the factor levels like this gets everything in the modeling to compare to control group!

NEP_range_model <- lmer(NEP_range ~ TREATMENT + (1|TANK_NUM), data=NEP_data2)
check_model(NEP_range_model)
summary(NEP_range_model)
anova(NEP_range_model)  # very sig effect p < 0.001 

emmeans(NEP_range_model, pairwise ~ "TREATMENT", adjust = "Tukey")

NEP_range_model_noRandom <- lm(NEP_range ~ TREATMENT, data = NEP_data2)
HSD.test(NEP_range_model_noRandom, "TREATMENT", console=TRUE)

# NEP mean #
NEP_plot <- chem_reframe_clean %>% # chem reframe data is only 12 and 9 pm sampling, CLEANED 
  ggplot(aes(x = TREATMENT, y = NEP_dailymean, color = TREATMENT)) +
  labs(x="",
       y = expression(bold("Daily Mean NEP" ~ (mmol ~ C ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_hline(yintercept = 0, lty=2) +
  geom_point(data = chem_reframe_clean, aes(x = TREATMENT, y = NEP_dailymean), alpha = 0.25) +
  stat_summary(size = 1, color = "black") + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1, color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")
NEP_plot
#ggsave(plot = NEP_plot, filename = here("Output", "NEP_Plots", "NEP_mean_plot.png"), width = 9, height = 9)

NEP_mean_model <- lmer(NEP_dailymean ~ TREATMENT + (1|TANK_NUM), data=chem_reframe_clean)
check_model(NEP_mean_model)
summary(NEP_mean_model)
anova(NEP_mean_model)

emmeans(NEP_mean_model, pairwise ~ "TREATMENT", adjust = "Tukey")

NEP_model_noRandom <- lm(NEP_dailymean ~ TREATMENT, data = chem_reframe_clean)
HSD.test(NEP_model_noRandom, "TREATMENT", console=TRUE)

####################################

## NEP DAY TIME ## 
NEP_data_day <- chem_reframe_DAY_clean %>% # plot NEP data during the DAY using cleaned data (cleaned 6/6/25)
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  select(DATE, TIME, TANK_NUM, TREATMENT, NEP) %>%
  reframe(NEP_day_mean = mean(NEP, na.rm = TRUE)) %>%
  drop_na()

NEP_day_treatment_means <- NEP_data_day %>%
  group_by(TREATMENT) %>%
  summarize(mean_NEP_treatment = mean(NEP_day_mean, na.rm = TRUE),
            se_NEP_treatment = sd(NEP_day_mean, na.rm = TRUE)/sqrt(n()))

NEP_data_day$TREATMENT <- factor(NEP_data_day$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEP_day_mean_plot <- NEP_data_day %>%
  ggplot(aes(x = TREATMENT, y = NEP_day_mean, color = TREATMENT)) +
  geom_hline(yintercept = 0, lty =2) +
  ylim(c(-3,5)) +
  labs(x = "", y = expression(bold("NCP" ~ (mmol ~ C ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = NEP_data_day, aes(x = TREATMENT, y = NEP_day_mean), alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  #geom_text(data = NEP_day_treatment_means, 
            #aes(x = TREATMENT, y = mean_NEP_treatment, 
               # label = paste0("", round(mean_NEP_treatment, 2), "±", round(se_NEP_treatment, 2))),
           # vjust = -7, hjust = 0.1, color = "black", size = 4) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none")
NEP_day_mean_plot
#ggsave(plot = NEP_day_mean_plot, filename = here("Output", "NEP_Plots", "NEP_day_mean_plot.png"), width = 9, height = 9)

NEP_data_day$TREATMENT <- factor(NEP_data_day$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEP_day_mean_model <- lmer(NEP_day_mean ~ TREATMENT + (1|TANK_NUM), data=NEP_data_day)
check_model(NEP_day_mean_model)
summary(NEP_day_mean_model)
anova(NEP_day_mean_model) # sig effect at p < 0.001

emmeans(NEP_day_mean_model, pairwise ~ "TREATMENT", adjust = "Tukey")
NEP_day_mean_model_noRandom <- lm(NEP_day_mean ~ TREATMENT, data = NEP_data_day)
HSD.test(NEP_day_mean_model_noRandom, "TREATMENT", console=TRUE)

## NEP NIGHT TIME ## 
NEP_data_night <- chem_reframe_NIGHT_clean %>% # plot NEP data during the NIGHT using cleaned data (cleaned 6/6/25)
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  select(DATE, TIME, TANK_NUM, TREATMENT, NEP) %>%
  reframe(NEP_night_mean = mean(NEP, na.rm = TRUE)) %>%
  drop_na()

NEP_night_treatment_means <- NEP_data_night %>%
  group_by(TREATMENT) %>%
  summarize(mean_NEP_treatment = mean(NEP_night_mean, na.rm = TRUE),
            se_NEP_treatment = sd(NEP_night_mean, na.rm = TRUE)/sqrt(n()))

NEP_data_night$TREATMENT <- factor(NEP_data_night$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEP_night_mean_plot <- NEP_data_night %>%
  ggplot(aes(x = TREATMENT, y = NEP_night_mean, color = TREATMENT)) +
  geom_hline(yintercept = 0, lty = 2) +
  ylim(c(-2,2)) +
  labs(x = "", y = "") +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = NEP_data_night, aes(x = TREATMENT, y = NEP_night_mean), alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", size = 3) + 
  #geom_text(data = NEP_night_treatment_means, 
           # aes(x = TREATMENT, y = mean_NEP_treatment, 
              #  label = paste0("", round(mean_NEP_treatment, 2), "±", round(se_NEP_treatment, 2))),
           # vjust = -7, hjust = 0.1, color = "black", size = 4) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none")
NEP_night_mean_plot
#ggsave(plot = NEP_night_mean_plot, filename = here("Output", "NEP_Plots", "NEP_night_mean_plot.png"), width = 9, height = 12)

NEP_data_night$TREATMENT <- factor(NEP_data_night$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEP_night_mean_model <- lmer(NEP_night_mean ~ TREATMENT + (1|TANK_NUM), data=NEP_data_night)
check_model(NEP_night_mean_model)
summary(NEP_night_mean_model)
anova(NEP_night_mean_model) # non sig effect of community on nighttime nep. p = 0.2065


## create patchwork plot
NEC_v_NEP_patchwork <- (NEC_day_mean_plot + NEC_night_mean_plot)/(NEP_day_mean_plot + NEP_night_mean_plot) + plot_annotation(tag_levels = "a")
NEC_v_NEP_patchwork
#ggsave(plot = NEC_v_NEP_patchwork, filename = here("Output", "NEC_v_NEP_patchwork.png"), width = 12, height = 10)


###########PH DATA ##########################
### PH DATA ###
pH_plot <- Data %>%
  group_by(TREATMENT, DATETIME)%>%
  summarise(mean_diff = mean(pHDiff, na.rm = TRUE),
            se_diff = sd(pHDiff, na.rm = TRUE)/sqrt(n()))%>%
  filter(!is.na(mean_diff), !is.na(se_diff), !is.na(DATETIME)) %>%
  ggplot(aes(x = DATETIME, y = mean_diff, color = TREATMENT, na.rm = TRUE))+
  geom_rect(data = rect_intervals, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 1/5, color = NA, inherit.aes = FALSE) + 
  scale_fill_identity() +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff+se_diff), color = "black", width = 0.1) + 
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  labs(x="Date & Time",
       y = "Change in pH") +
  geom_line() +
  annotate("text", x = ymd_hms("2024-06-03 08:00:00"), y = 0.1, label = "Overcast", size = 5) +
  annotate("text", x = ymd_hms("2024-06-06 13:00:00"), y = 0.1, label = "Overcast \n Rain", size = 5) +
  theme(plot.title = element_text(size = 14))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
pH_plot

# create columns for mean DAY pH and mean NIGHT pH # 
pH_day_night <- Data %>% 
  filter(TIME %in% c("12:00:00", "21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(pH_day = pH[TIME == hms("12:00:00")], 
          pH_night = pH[TIME == hms("21:00:00")])
#write_csv(pH_day_night, here("Data", "Chemistry", "pH_day_night.csv"))

# now calculate mean pH for both day and night sampling times # 
pH_day_night_means <- pH_day_night %>% 
  group_by(TREATMENT) %>%
  summarize(pH_day_mean = mean(pH_day, na.rm = TRUE), 
            pH_day_se = sd(pH_day, na.rm = TRUE)/sqrt(n()), 
            pH_night_mean = mean(pH_night, na.rm = TRUE), 
            pH_night_se = sd(pH_night, na.rm = TRUE)/sqrt(n()))
#write_csv(pH_day_night_means, here("Data", "Chemistry", "pH_day_night_means.csv"))

# pH daytime mean #
pH_day_night_means$TREATMENT <- factor(pH_day_night_means$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))
pHday_means_plot <- ggplot(data = pH_day_night_means, aes(x = TREATMENT, y = pH_day_mean, color = TREATMENT)) +
  geom_point(size = 2.5) + 
  labs(x = "Community Tank", y = "Daytime Mean pH") + 
  geom_jitter(data = pH_day_night, aes(x = TREATMENT, y = pH_day), alpha = 0.4) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble/CCA-Dominated")) +
  theme(plot.title = element_text(size = 14))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
pHday_means_plot
#ggsave(plot = pHday_means_plot, filename = here("Output", "pHPlots", "pHday_means_plot.png"), width = 9, height = 9)


# pH nighttime mean # 
pHnight_means_plot <- ggplot(data = pH_day_night_means, aes(x = TREATMENT, y = pH_night_mean, color = TREATMENT)) +
  geom_point(size = 2.5) + 
  labs(x = "Community Tank", y = "Nighttime Mean pH") + 
  geom_jitter(data = pH_day_night, aes(x = TREATMENT, y = pH_night), alpha = 0.4) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble/CCA-Dominated")) +
  theme(plot.title = element_text(size = 14))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"),
                     values = c("blue", "darkgreen", "coral", "tan"))
pHnight_means_plot
#ggsave(plot = pHnight_means_plot, filename = here("Output", "pHPlots", "pHnight_means_plot.png"), width = 9, height = 9)


## filter only 12:00 and 21:00 sampling times ## 
pH_clean <- Data %>% 
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(pH_range = pH[TIME == hms("12:00:00")] - pH[TIME == hms("21:00:00")],
            pH_dailymean = mean(pH, na.rm = TRUE))
pH_clean

# add pH day and pH night means to pH_clean 
pH_clean <- pH_clean %>% 
  right_join(pH_day_night_means)

#write_csv(pH_clean, here("Data", "Chemistry", "Cleaned_pH_Data_FULL.csv"))

pH_plotdata<- pH_clean %>%
  group_by(TREATMENT) %>%
  summarize(pH_rangemean = mean(pH_range, na.rm = TRUE),
            pH_rangese = sd(pH_range, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH_dailymean, na.rm = TRUE),
            pH_se = sd(pH_dailymean, na.rm = TRUE)/sqrt(n()))
pH_plotdata
#write_csv(pH_plotdata, here("Data", "Chemistry", "Cleaned_pH_Data_per_Treatment.csv"))

pH_plotdata_tank<- pH_clean %>%
  group_by(TANK_NUM) %>%
  summarize(pH_rangemean = mean(pH_range, na.rm = TRUE),
            pH_rangese = sd(pH_range, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH_dailymean, na.rm = TRUE),
            pH_se = sd(pH_dailymean, na.rm = TRUE)/sqrt(n()))
pH_plotdata_tank

#write_csv(pH_plotdata_tank_full, here("Data", "Chemistry", "Cleaned_pH_Data_per_Tank.csv"))


## plot pH range from 12:00 and 21:00 sampling throughout experiment ## 
pH_range_plot <- chem_reframe_clean %>%
  ggplot(aes(x = TREATMENT, y = pH_range, color = TREATMENT)) +
  geom_point(data = chem_reframe_clean, aes(x = TREATMENT, y = pH_range), alpha = 0.25) +
  labs(x = "Treatment", y = "Daily Range pH") +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1, color = "black") +
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") 
pH_range_plot


pH_range_model <- lmer(pH_range ~ TREATMENT +(1|TANK_NUM), data= chem_reframe_clean)
check_model(pH_range_model)
summary(pH_range_model)
anova(pH_range_model) # significant effect of treatment on the range in pH. all treatments sig diff from control

emmeans(pH_range_model, pairwise ~ "TREATMENT", adjust = "Tukey")
pH_range_model_noRandom <- lm(pH_range ~ TREATMENT, data = chem_reframe_clean)
HSD.test(pH_range_model_noRandom, "TREATMENT", console=TRUE)

## mean pH plot of 12:00 and 21:00 sampling throughout experiment ## 
chem_reframe_clean$TREATMENT <- factor(chem_reframe_clean$TREATMENT, levels = c("Control", "Algae_Dom", 
                                                                                "Coral_Dom", "Rubble_Dom"))
pH_mean_plot <- chem_reframe_clean %>%
  ggplot(aes(x = TREATMENT, y = pH_dailymean, color = TREATMENT)) +
  geom_point(data = chem_reframe_clean, aes(x = TREATMENT, y = pH_dailymean), alpha = 0.25) +
  labs(x = "Treatment", y = "Daily Mean pH") +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1, color = "black") +
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") 
pH_mean_plot
#ggsave(plot = pH_mean_plot, filename = here("Output", "pHPlots", "pH_mean_plot.png"), width = 9, height = 6)

## mean pH treatment stats ## 
mean_pH_model <- lmer(pH_dailymean ~ TREATMENT +(1|TANK_NUM), data= chem_reframe_clean)
check_model(mean_pH_model)
summary(mean_pH_model)
anova(mean_pH_model)
# no significant effect of community type on daily mean pH 

# mean DOC and treatments #
DOC_mean_plot <- chem_reframe_clean %>%
  ggplot(aes(x = TREATMENT, y = DOC_dailymean, color = TREATMENT)) +
  geom_point(data = chem_reframe_clean, aes(x = TREATMENT, y = DOC_dailymean), alpha = 0.25) +
  labs(x = "Treatment", y = "Daily Mean DOC") +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1, color = "black") +
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") 
DOC_mean_plot

mean_DOC_model <- lmer(DOC_dailymean ~ TREATMENT +(1|TANK_NUM), data= chem_reframe_clean)
check_model(mean_DOC_model)
summary(mean_DOC_model)
anova(mean_DOC_model)

# range DOC and treatment 
DOC_range_plot <- chem_reframe_clean %>%
  ggplot(aes(x = TREATMENT, y = DOC_range, color = TREATMENT)) +
  geom_point(data = chem_reframe_clean, aes(x = TREATMENT, y = DOC_range), alpha = 0.25) +
  labs(x = "Treatment", y = "Daily Range DOC") +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1, color = "black") +
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") 
DOC_range_plot

range_DOC_model <- lmer(DOC_range ~ TREATMENT +(1|TANK_NUM), data = chem_reframe_clean)
check_model(range_DOC_model)
summary(range_DOC_model)
anova(range_DOC_model)



############# NEP DATA ############

#### Look at NEP vs pH and DOC
Clean_Chem_all<- chem_reframe_NIGHT_clean %>% 
  bind_rows(chem_reframe_DAY_clean) %>%
  filter(NEP < 4) # remove the one outlier with really high NEP at night

# ANCOVA with pH~NEP*Treament and accounting for tankID
mod_pH<-lm(pH ~ NEP * TREATMENT, data = Clean_Chem_all)
anova(mod_pH)
summary(mod_pH)

mod_DOC<-lm(DOC~NEP*TREATMENT,data = Clean_Chem_all)
anova(mod_DOC)
summary(mod_DOC)

mod_DOC_NEC<-lm(DOC~NEC*TREATMENT,data = Clean_Chem_all)
anova(mod_DOC_NEC)
summary(mod_DOC_NEC)

mod_NEP_NEC<-lm(NEC~pH*TREATMENT,data = Clean_Chem_all)
anova(mod_NEP_NEC)
summary(mod_NEP_NEC)

chem_plot_all_pH<-Clean_Chem_all %>%
  ggplot(aes(x = NEP, y = pH))+
  geom_point(aes(color = TREATMENT), alpha = 0.2)+
  geom_smooth(method = "lm", color = "black")+
  labs(color = "",
       y = expression("pH"[T]),
       x = "")+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated",
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"), guide = FALSE)
chem_plot_all_pH

chem_plot_all_DOC <- Clean_Chem_all %>%
  ggplot(aes(x = NEP, y = DOC))+
  geom_point(aes(color = TREATMENT), alpha = 0.2)+
  geom_smooth(method = "lm", color = "black")+
  labs(color = "",
       x = expression("NEP (mmol C m"^2~" hr"^-1~")"),
       y = expression("DOC ("~mu~"mol L"^-1~")"))+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated",
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"))

chem_plot_all_DOC

chem_plot_all_NEP_NEC<-Clean_Chem_all %>%
  ggplot(aes(x = NEP, y = NEC, color = TREATMENT))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "lm")+
  labs(color = "",
       y = expression("NEC (mmol C m"^2~" hr"^-1~")"),
       x = "")+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        # legend.position = "bottom",
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated",
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"), guide = FALSE)

chem_plot_all_NEP_NEC

(chem_plot_all_NEP_NEC/chem_plot_all_pH/chem_plot_all_DOC& theme(legend.position = 'bottom'))+
  plot_annotation(tag_levels = "a")+plot_layout(guides = "collect")


##### NEP and pH #####
Data1 <- Data %>%
  filter(!NEP < -4) %>%
  filter(!TA > 2700) %>%
  filter(!pH < 7.8)

chem_reframe_clean1$TREATMENT <- factor(chem_reframe_clean1$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))
chem_reframe_clean <- chem_reframe_clean %>%
  filter(!deltapH_dailymean > 0.2,
         !deltapH_dailymean < -0.1,
         !pH_dailymean < 7.99,
         !pH_dailymean > 8.2)

NEP_deltapH_plot <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEP_dailymean, y = deltapH_dailymean)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1, label.y = 0.15, size = 5) + 
  stat_cor(label.x = 1, label.y = 0.125, size = 5)
NEP_deltapH_plot
#ggsave(plot = NEP_deltapH_plot, filename = here("Output", "NEP_Plots", "NEP_deltapH.png"), width = 12, height = 9)

NEP_deltapH_model <- lmer(deltapH_dailymean ~ NEP_dailymean * TREATMENT + (1|TANK_NUM), data = chem_reframe_clean1)
check_model(NEP_deltapH_model)
summary(NEP_deltapH_model)
###########################################
NEP_deltapH_plot2 <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEP_dailymean, y = deltapH_dailymean)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "NEP", y = "deltapH") +
  facet_wrap(~TREATMENT) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1, label.y = 0.15, size = 5) + 
  stat_cor(label.x = 1, label.y = 0.125, size = 5)
NEP_deltapH_plot2
#ggsave(plot = NEP_deltapH_plot2, filename = here("Output", "NEP_Plots", "NEP_deltapH_facet.png"), width = 12, height = 9)
###########################################

##### NEP and DIC ##### 

NEP_DIC_plot <- Data1 %>%
  ggplot(aes(x = NEP, y = DIC_µmol_kg)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 3, label.y = 2100, size = 5) + 
  stat_cor(label.x = 3, label.y = 2070, size = 5)
NEP_DIC_plot
#ggsave(plot = NEP_DIC_plot, filename = here("Output", "NEP_Plots", "NEP_DIC.png"), width = 12, height = 9)


##### NEC and pH ##### 

NEC_deltapH_plot <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEC_dailymean, y = deltapH_dailymean)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "∆ pH") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = -1, label.y = 0.15, size = 5) + 
  stat_cor(label.x = -1, label.y = 0.125, size = 5) 
NEC_deltapH_plot
#ggsave(plot = NEC_deltapH_plot, filename = here("Output", "TA_NECPlots", "NEC_deltapH.png"), width = 12, height = 9)

NEC_deltapH_model <- lmer(deltapH_dailymean ~ NEC_dailymean * TREATMENT + (1|TANK_NUM), data = chem_reframe_clean1)
check_model(NEC_deltapH_model)
summary(NEC_deltapH_model)

######################################################
NEC_deltapH_plot2 <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEC_dailymean, y = deltapH_dailymean)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "NEC", y = "∆ pH") +
  facet_wrap(~TREATMENT) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = -1, label.y = 0.15, size = 5) + 
  stat_cor(label.x = -1, label.y = 0.125, size = 5) 
NEC_deltapH_plot2
#ggsave(plot = NEC_deltapH_plot2, filename = here("Output", "TA_NECPlots", "NEC_deltapH_facet.png"), width = 12, height = 9)

#################################################################
##### NEC and DIC ##### 
NEC_DIC_plot <- Data1 %>%
  ggplot(aes(x = NEC, y = DIC_µmol_kg)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "DIC (µmol/kg)") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1, label.y = 2100, size = 5) + 
  stat_cor(label.x = 1, label.y = 2070, size = 5)
NEC_DIC_plot
#ggsave(plot = NEC_DIC_plot, filename = here("Output", "TA_NECPlots", "NEC_DIC.png"), width = 12, height = 9)


#######################################
##### NEC AND TA #####
NEC_TA_plot <- Data1 %>%
  ggplot(aes(x = NEC, y = TA)) +
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "NEC", y = "TA") +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 2.5, label.y = 2390, size = 5) + 
  stat_cor(label.x = 2.5, label.y = 2375, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"))
NEC_TA_plot
#ggsave(plot = NEC_TA_plot, filename = here("Output", "TA_NECPlots", "NEC_TA.png"), width = 12, height = 9)

NEC_TA_model <- lmer(TA ~ NEC + (1|TANK_NUM), data = Data1)
check_model(NEC_TA_model)
summary(NEC_TA_model)

NEC_TA_facetplot <- Data1 %>%
  ggplot(aes(x = NEC, y = TA)) +
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) +
  facet_wrap(~ TREATMENT) +
  labs(x = "NEC", y = "TA") +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 1, label.y = 2420, size = 5) + 
  stat_cor(label.x = 1, label.y = 2390, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"))
NEC_TA_facetplot
#ggsave(plot = NEC_TA_facetplot, filename = here("Output", "TA_NECPlots", "NEC_TA_treatments.png"), width = 14, height = 9)

##### NEC and DOC ##### 
NEC_deltaDOC_plot <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEC_dailymean, y = log(deltaDOC_dailymean))) + # log scaled to improve normality
  geom_point(aes(color = TREATMENT)) + 
  xlim(c(-1,4)) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1))), y = expression(bold("∆ DOC"))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1.5, label.y = 0.75, size = 5) + 
  stat_cor(label.x = 1.5, label.y = 0, size = 5)
NEC_deltaDOC_plot
#ggsave(plot = NEC_deltaDOC_plot, filename = here("Output", "TA_NECPlots", "NEC_deltaDOC.png"), width = 14, height = 10)

NEC_deltaDOC_model <- lmer(log(deltaDOC_dailymean) ~ NEC_dailymean + (1|TANK_NUM), data = chem_reframe_clean1)
check_model(NEC_deltaDOC_model)
summary(NEC_deltaDOC_model) 

#################################################
NEC_deltaDOC_plot2 <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEC_dailymean, y = log(deltaDOC_dailymean))) + # log scaled to improve normality
  geom_point(aes(color = TREATMENT)) + 
  xlim(c(-1,4)) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1))), y = expression(bold("∆ DOC"))) +
  facet_wrap(~TREATMENT) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1.5, label.y = 0.75, size = 5) + 
  stat_cor(label.x = 1.5, label.y = 0, size = 5)
NEC_deltaDOC_plot2
#ggsave(plot = NEC_deltaDOC_plot2, filename = here("Output", "TA_NECPlots", "NEC_deltaDOC_facet.png"), width = 14, height = 10)
#######################################################


##### NEP and DOC ##### 
NEP_deltaDOC_plot <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEP_dailymean, y = log(deltaDOC_dailymean))) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEP" ~ (mmol ~ C ~ m^2 ~ h^-1))), y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1.5, label.y = 1, size = 5) + 
  stat_cor(label.x = 1.5, label.y = 0.55, size = 5)
NEP_deltaDOC_plot
#ggsave(plot = NEP_deltaDOC_plot, filename = here("Output", "NEP_Plots", "NEP_deltaDOC.png"), width = 14, height = 10)

NEP_deltaDOC_model <- lmer(log(deltaDOC_dailymean) ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean)
check_model(NEP_deltaDOC_model)
summary(NEP_deltaDOC_model)

##############################################
NEP_deltaDOC_plot2 <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEP_dailymean, y = log(deltaDOC_dailymean))) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEP" ~ (mmol ~ C ~ m^2 ~ h^-1))), y = "") +
  facet_wrap(~TREATMENT) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1.5, label.y = 1, size = 5) + 
  stat_cor(label.x = 1.5, label.y = 0.55, size = 5)
NEP_deltaDOC_plot2
#ggsave(plot = NEP_deltaDOC_plot2, filename = here("Output", "NEP_Plots", "NEP_deltaDOC_facet.png"), width = 14, height = 10)
###########################################

NEC_NEP_biogeochem <- (NEC_deltapH_plot + NEP_deltapH_plot)/(NEC_deltaDOC_plot + NEP_deltaDOC_plot) + plot_annotation(tag_levels = "a")
NEC_NEP_biogeochem
#ggsave(plot = NEC_NEP_biogeochem, filename = here("Output", "NEC_NEP_biogeochem_patch.png"), width = 15, height = 15)

### NEP and NEC vs DIC and DOC ### 
NEC_NEP_deltapH_DOC <- (NEC_deltapH_plot + NEP_deltapH_plot)/(NEC_deltaDOC_plot + NEP_deltaDOC_plot) + plot_annotation(tag_levels = "a")
NEC_NEP_deltapH_DOC
#ggsave(plot = NEC_NEP_deltapH_DOC, filename = here("Output", "NEC_NEP_deltapH_DOC_patch.png"), width = 15, height = 15)

############### nep and nec with raw data ############
# NEP and pH plot - raw data # 
chem_reframe_clean1 <- chem_reframe_clean %>%
  filter(!pH_dailymean > 8.2,
         !NEP_dailymean > 4.0, 
         !NEP_dailymean < -2)


NEP_pH_plot <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEP_dailymean, y = pH_dailymean)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1.5, label.y = 8.12, size = 5) + 
  stat_cor(label.x = 1.5, label.y = 8.10, size = 5)
NEP_pH_plot
#ggsave(plot = NEP_pH_plot, filename = here("Output", "NEP_Plots", "NEP_pH.png"), width = 12, height = 9)

NEP_pH_plot_model <- lmer(pH_dailymean ~ NEP_dailymean * TREATMENT + (1|TANK_NUM), data = chem_reframe_clean)
check_model(NEP_pH_plot_model)
plot(fitted(NEP_pH_plot_model), residuals(NEP_pH_plot_model),
     xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
summary(NEP_pH_plot_model)

###########################################
NEP_pH_plot2 <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEP_dailymean, y = pH_dailymean)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "") +
  facet_wrap(~TREATMENT) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1.5, label.y = 8.12, size = 5) + 
  stat_cor(label.x = 1.5, label.y = 8.10, size = 5)
NEP_pH_plot2
#ggsave(plot = NEP_pH_plot2, filename = here("Output", "NEP_pH_facet.png"), width = 12, height = 9)
#########################################

# NEC and pH plot - raw data #
NEC_pH_plot <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEC_dailymean, y = pH_dailymean)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "Mean pH") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = -3, label.y = 8.10, size = 5) + 
  stat_cor(label.x = -3, label.y = 8.08, size = 5) 
NEC_pH_plot
#ggsave(plot = NEC_pH_plot, filename = here("Output", "TA_NECPlots", "NEC_pH.png"), width = 12, height = 9)

NEC_pH_model <- lmer(pH_dailymean ~ NEC_dailymean * TREATMENT + (1|TANK_NUM), data = chem_reframe_clean)
check_model(NEC_pH_model)
plot(fitted(NEC_pH_model), residuals(NEC_pH_model),
     xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
summary(NEC_pH_model) 

############################################
NEC_pH_plot2 <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEC_dailymean, y = pH_dailymean)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "Mean pH") +
  facet_wrap(~TREATMENT) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = -3, label.y = 8.10, size = 5) + 
  stat_cor(label.x = -3, label.y = 8.08, size = 5) 
NEC_pH_plot2
#ggsave(plot = NEC_pH_plot2, filename = here("Output", "NEC_pH_facet.png"), width = 12, height = 9)
############################################

# NEC and DOC plot - raw data #
NEC_DOC_plot <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEC_dailymean, y = log(DOC_dailymean))) + # log scaled to improve normality
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NCC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1))), y = expression(bold("Mean DOC" ~ (log[10])))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = -3, label.y = 5.5, size = 5) + 
  stat_cor(label.x = -3, label.y = 5.3, size = 5)
NEC_DOC_plot
#ggsave(plot = NEC_DOC_plot, filename = here("Output", "TA_NECPlots", "NEC_DOC.png"), width = 14, height = 10)

NEC_DOC_model <- lmer(log(DOC_dailymean) ~ NEC_dailymean * TREATMENT + (1|TANK_NUM), data = chem_reframe_clean)
check_model(NEC_DOC_model)
summary(NEC_DOC_model) # p = 0.096

##########################################
NEC_DOC_plot2 <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEC_dailymean, y = log(DOC_dailymean))) + # log scaled to improve normality
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NCC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1))), y = expression(bold("Mean DOC" ~ (log[10])))) +
  facet_wrap(~TREATMENT) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = -3, label.y = 5.5, size = 5) + 
  stat_cor(label.x = -3, label.y = 5.3, size = 5)
NEC_DOC_plot2
#ggsave(plot = NEC_DOC_plot2, filename = here("Output", "NEC_DOC_facet.png"), width = 14, height = 10)
##########################################

# NEP and DOC plot - raw data # 
NEP_DOC_plot <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEP_dailymean, y = log(DOC_dailymean))) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NCP" ~ (mmol ~ C ~ m^2 ~ h^-1))), y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1, label.y = 5.8, size = 5) + 
  stat_cor(label.x = 1, label.y = 5.6, size = 5)
NEP_DOC_plot
#ggsave(plot = NEP_DOC_plot, filename = here("Output", "NEP_Plots", "NEP_DOC.png"), width = 14, height = 10)

NEP_DOC_model <- lmer(log(DOC_dailymean) ~ NEP_dailymean * TREATMENT + (1|TANK_NUM), data = chem_reframe_clean)
check_model(NEP_DOC_model)
summary(NEP_DOC_model) # sig effect of NEP on mean DOC. p = 0.00684.

##########################################
NEP_DOC_plot2 <- chem_reframe_clean1 %>%
  ggplot(aes(x = NEP_dailymean, y = log(DOC_dailymean))) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NCP" ~ (mmol ~ C ~ m^2 ~ h^-1))), y = "") +
  facet_wrap(~TREATMENT) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1, label.y = 5.8, size = 5) + 
  stat_cor(label.x = 1, label.y = 5.6, size = 5)
NEP_DOC_plot2
#ggsave(plot = NEP_DOC_plot2, filename = here("Output", "NEP_DOC_facet.png"), width = 14, height = 10)


# patchwork plot of daily mean NEC and NEP on daily mean pH and DOC #
NEC_NEP_raw_biogeochem <- (NEC_pH_plot + NEP_pH_plot)/(NEC_DOC_plot + NEP_DOC_plot) + plot_annotation(tag_levels = "a")
NEC_NEP_raw_biogeochem
#ggsave(plot = NEC_NEP_raw_biogeochem, filename = here("Output", "NEC_NEP_raw_biogeochem_patch.png"), width = 15, height = 15)

## NEP/NEC and biogeochem ancova ## 
chem_reframe_clean$TREATMENT <- as.factor(chem_reframe_clean$TREATMENT)

pH_nep_aov <- aov(NEP_dailymean ~ TREATMENT, data = chem_reframe_clean)
summary(pH_nep_aov) # p < 0.05, means that the covariate (NEP daily mean) and TREATMENT are NOT independent
leveneTest(pH_dailymean ~ TREATMENT, data = chem_reframe_clean) # p > 0.05 means that the variances between groups are equal (homogeneity of variance)

pH_nep_ancova <- lm(pH_dailymean ~ TREATMENT + NEP_dailymean, data = chem_reframe_clean)
check_model(pH_nep_ancova)
plot(pH_dailymean ~ NEP_dailymean, data = chem_reframe_clean)
anova(pH_nep_ancova)

pH_nec_aov <- aov(NEC_dailymean ~ TREATMENT, data = chem_reframe_clean)
summary(pH_nec_aov) # p < 0.05, means that the covariate (NEC daily mean) and TREATMENT are NOT independent

pH_nec_ancova <- lm(pH_dailymean ~ TREATMENT + NEC_dailymean, data = chem_reframe_clean)
check_model(pH_nec_ancova)
plot(pH_dailymean ~ NEC_dailymean, data = chem_reframe_clean)
anova(pH_nec_ancova)

leveneTest(DOC_dailymean ~ TREATMENT, data = chem_reframe_clean) # homogeneity of variances 

doc_nep_ancova <- lm(log(DOC_dailymean) ~ TREATMENT + NEP_dailymean, data = chem_reframe_clean)
check_model(doc_nep_ancova)
plot(log(DOC_dailymean) ~ NEP_dailymean, data = chem_reframe_clean) 
anova(doc_nep_ancova) ## sig effect of TREATMENT p = 0.036149, and NEP daily mean p = 0.004326. 

doc_nec_ancova <- lm(log(DOC_dailymean) ~ TREATMENT + NEC_dailymean, data = chem_reframe_clean)
check_model(doc_nec_ancova)
plot(log(DOC_dailymean) ~ NEC_dailymean, data = chem_reframe_clean) 
anova(doc_nec_ancova) ## sig effect of TREATMENT p = 0.04775


##### NEC and NEP v delta pH and delta DOC PER TREATMENT #####
### ALGAE DOM ###
# NEP v delta pH #
chem_reframe_clean_ALGAE <- chem_reframe_clean1 %>%
  filter(TREATMENT == "Algae_Dom") %>%
  select(TANK_NUM, TA_range:NEP_dailymean)
chem_reframe_clean_ALGAE

chem_reframe_clean_CORAL <- chem_reframe_clean1 %>%
  filter(TREATMENT == "Coral_Dom") %>%
  select(TANK_NUM, TA_range:NEP_dailymean)
chem_reframe_clean_CORAL

chem_reframe_clean_RUBBLECCA <- chem_reframe_clean1 %>%
  filter(TREATMENT == "Rubble_Dom") %>%
  select(TANK_NUM, TA_range:NEP_dailymean)
chem_reframe_clean_RUBBLECCA

chem_reframe_clean_CONTROL <- chem_reframe_clean1 %>%
  filter(TREATMENT == "Control") %>%
  select(TANK_NUM, TA_range:NEP_dailymean)
chem_reframe_clean_CONTROL

NEP_deltapH_algae <- chem_reframe_clean_ALGAE %>%
  ggplot(aes(x = NEP_dailymean, y = deltapH_dailymean)) +
  geom_point(color = "darkgreen") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1, label.y = -0.01, size = 5) + 
  stat_cor(label.x = 1, label.y = -0.02, size = 5)
NEP_deltapH_algae

NEP_deltapH_algae_model <- lmer(deltapH_dailymean ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean_ALGAE)
check_model(NEP_deltapH_algae_model)
summary(NEP_deltapH_algae_model) # p = 0.05

# NEP v delta DOC #
NEP_deltaDOC_algae <- chem_reframe_clean_ALGAE %>%
  ggplot(aes(x = NEP_dailymean, y = log(deltaDOC_dailymean))) +
  geom_point(color = "darkgreen") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEP" ~ (mmol ~ C ~ m^2 ~ h^-1))), y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1.5, label.y = 2, size = 5) + 
  stat_cor(label.x = 1.5, label.y = 1.5, size = 5)
NEP_deltaDOC_algae

NEP_deltaDOC_algae_model <- lmer(log(deltaDOC_dailymean) ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean_ALGAE)
check_model(NEP_deltaDOC_algae_model)
summary(NEP_deltaDOC_algae_model)

# NEC v delta pH # 
NEC_deltapH_algae <- chem_reframe_clean_ALGAE %>%
  ggplot(aes(x = NEC_dailymean, y = deltapH_dailymean)) +
  geom_point(color = "darkgreen") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "∆ pH") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = -1.5, label.y = -0.025, size = 5) + 
  stat_cor(label.x = -1.5, label.y = -0.030, size = 5)
NEC_deltapH_algae

NEC_deltapH_algae_model <- lmer(deltapH_dailymean ~ NEC_dailymean + (1|TANK_NUM), data = chem_reframe_clean_ALGAE)
check_model(NEC_deltapH_algae_model)
summary(NEC_deltapH_algae_model)

# NEC v delta DOC # 
NEC_deltaDOC_algae <- chem_reframe_clean_ALGAE %>%
  ggplot(aes(x = NEC_dailymean, y = log(deltaDOC_dailymean))) +
  geom_point(color = "darkgreen") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1))), y = expression(bold("∆ DOC"))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = -1.5, label.y = 2, size = 5) + 
  stat_cor(label.x = -1.5, label.y = 1.5, size = 5)
NEC_deltaDOC_algae

NEC_deltaDOC_algae_model <- lmer(log(deltaDOC_dailymean) ~ NEC_dailymean + (1|TANK_NUM), data = chem_reframe_clean_ALGAE)
check_model(NEC_deltaDOC_algae_model)
summary(NEC_deltaDOC_algae_model)

## plot for ALGAE ##
NEC_NEP_pH_DOC_ALGAE <- (NEC_deltapH_algae + NEP_deltapH_algae)/(NEC_deltaDOC_algae + NEP_deltaDOC_algae) + plot_annotation(tag_levels = "a")
NEC_NEP_pH_DOC_ALGAE
#ggsave(plot = NEC_NEP_pH_DOC_ALGAE, filename = here("Output", "NEC_NEP_pH_DOC_ALGAEpatch.png"), width = 15, height = 15)

### CORAL DOM ###
# NEP v delta pH #
NEP_deltapH_coral <- chem_reframe_clean_CORAL %>%
  ggplot(aes(x = NEP_dailymean, y = deltapH_dailymean)) +
  geom_point(color = "coral") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
 stat_regline_equation(label.x = -0.5, label.y = -0.02, size = 5) + 
 stat_cor(label.x = -0.5, label.y = -0.03, size = 5)
NEP_deltapH_coral

NEP_deltapH_coral_model <- lmer(deltapH_dailymean ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean_CORAL)
check_model(NEP_deltapH_coral_model)
summary(NEP_deltapH_coral_model)

# NEP v delta DOC #
NEP_deltaDOC_coral <- chem_reframe_clean_CORAL %>%
  ggplot(aes(x = NEP_dailymean, y = log(deltaDOC_dailymean))) +
  geom_point(color = "coral") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEP" ~ (mmol ~ C ~ m^2 ~ h^-1))), y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
 stat_regline_equation(label.x = -1, label.y = 0.5, size = 5) + 
 stat_cor(label.x = -1, label.y = 0.25, size = 5)
NEP_deltaDOC_coral

NEP_deltaDOC_coral_model <- lmer(log(deltaDOC_dailymean) ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean_CORAL)
check_model(NEP_deltaDOC_coral_model)
summary(NEP_deltaDOC_coral_model)

# NEC v delta pH # 
NEC_deltapH_coral <- chem_reframe_clean_CORAL %>%
  ggplot(aes(x = NEC_dailymean, y = deltapH_dailymean)) +
  geom_point(color = "coral") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "∆ pH") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
 stat_regline_equation(label.x = -1, label.y = -0.02, size = 5) + 
 stat_cor(label.x = -1, label.y = -0.03, size = 5)
NEC_deltapH_coral

NEC_deltapH_coral_model <- lmer(deltapH_dailymean ~ NEC_dailymean + (1|TANK_NUM), data = chem_reframe_clean_CORAL)
check_model(NEC_deltapH_coral_model)
summary(NEC_deltapH_coral_model)

# NEC v delta DOC # 
NEC_deltaDOC_coral <- chem_reframe_clean_CORAL %>%
  ggplot(aes(x = NEC_dailymean, y = log(deltaDOC_dailymean))) +
  geom_point(color = "coral") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1))), y = expression(bold("∆ DOC"))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
 stat_regline_equation(label.x = -1, label.y = 1, size = 5) + 
 stat_cor(label.x = -1, label.y = 0.5, size = 5)
NEC_deltaDOC_coral

NEC_deltaDOC_coral_model <- lmer(log(deltaDOC_dailymean) ~ NEC_dailymean + (1|TANK_NUM), data = chem_reframe_clean_CORAL)
check_model(NEC_deltaDOC_coral_model)
summary(NEC_deltaDOC_coral_model)

## plot for CORAL ##
NEC_NEP_pH_DOC_CORAL <- (NEC_deltapH_coral + NEP_deltapH_coral)/(NEC_deltaDOC_coral + NEP_deltaDOC_coral) + plot_annotation(tag_levels = "a")
NEC_NEP_pH_DOC_CORAL
#ggsave(plot = NEC_NEP_pH_DOC_CORAL, filename = here("Output", "NEC_NEP_pH_DOC_CORALpatch.png"), width = 15, height = 15)

### CONTROL ###
# NEP v delta pH #
NEP_deltapH_control <- chem_reframe_clean_CONTROL %>%
  ggplot(aes(x = NEP_dailymean, y = deltapH_dailymean)) +
  geom_point(color = "blue") + 
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
 stat_regline_equation(label.x = 0.5, label.y = 0, size = 5) + 
 stat_cor(label.x = 0.5, label.y = -0.003, size = 5)
NEP_deltapH_control

NEP_deltapH_control_model <- lmer(deltapH_dailymean ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean_CONTROL)
check_model(NEP_deltapH_control_model)
summary(NEP_deltapH_control_model) # p = 0.0017

# NEP v delta DOC #
NEP_deltaDOC_control <- chem_reframe_clean_CONTROL %>%
  ggplot(aes(x = NEP_dailymean, y = log(deltaDOC_dailymean))) +
  geom_point(color = "blue") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEP" ~ (mmol ~ C ~ m^2 ~ h^-1))), y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
 stat_regline_equation(label.x = 0.5, label.y = 3.5, size = 5) + 
 stat_cor(label.x = 0.5, label.y = 3.25, size = 5)
NEP_deltaDOC_control

NEP_deltaDOC_control_model <- lmer(log(deltaDOC_dailymean) ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean_CONTROL)
check_model(NEP_deltaDOC_control_model)
summary(NEP_deltaDOC_control_model)

# NEC v delta pH # 
NEC_deltapH_control <- chem_reframe_clean_CONTROL %>%
  ggplot(aes(x = NEC_dailymean, y = deltapH_dailymean)) +
  geom_point(color = "blue") + 
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "∆ pH") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
 stat_regline_equation(label.x = -0.8, label.y = 0.005, size = 5) + 
 stat_cor(label.x = -0.8, label.y = 0, size = 5)
NEC_deltapH_control

NEC_deltapH_control_model <- lmer(deltapH_dailymean ~ NEC_dailymean + (1|TANK_NUM), data = chem_reframe_clean_CONTROL)
check_model(NEC_deltapH_control_model)
summary(NEC_deltapH_control_model) # p = 0.000116

# NEC v delta DOC # 
NEC_deltaDOC_control <- chem_reframe_clean_CONTROL %>%
  ggplot(aes(x = NEC_dailymean, y = log(deltaDOC_dailymean))) +
  geom_point(color = "blue") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1))), y = expression(bold("∆ DOC"))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
 stat_regline_equation(label.x = -0.5, label.y = 3.5, size = 5) + 
 stat_cor(label.x = -0.5, label.y = 3.25, size = 5)
NEC_deltaDOC_control

NEC_deltaDOC_control_model <- lmer(log(deltaDOC_dailymean) ~ NEC_dailymean + (1|TANK_NUM), data = chem_reframe_clean_CONTROL)
check_model(NEC_deltaDOC_control_model)
summary(NEC_deltaDOC_control_model)

## plot for control ##
NEC_NEP_pH_DOC_CONTROL <- (NEC_deltapH_control + NEP_deltapH_control)/(NEC_deltaDOC_control + NEP_deltaDOC_control) + plot_annotation(tag_levels = "a")
NEC_NEP_pH_DOC_CONTROL
#ggsave(plot = NEC_NEP_pH_DOC_CONTROL, filename = here("Output", "NEC_NEP_pH_DOC_CONTROLpatch.png"), width = 15, height = 15)

### RUBBLE/CCA DOM ###
# NEP v delta pH #
NEP_deltapH_rubble <- chem_reframe_clean_RUBBLECCA %>%
  ggplot(aes(x = NEP_dailymean, y = deltapH_dailymean)) +
  geom_point(color = "tan") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1, label.y = 0.05, size = 5) + 
  stat_cor(label.x = 1, label.y = 0.045, size = 5)
NEP_deltapH_rubble

NEP_deltapH_rubble_model <- lmer(deltapH_dailymean ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean_RUBBLECCA)
check_model(NEP_deltapH_rubble_model)
summary(NEP_deltapH_rubble_model)

# NEP v delta DOC #
NEP_deltaDOC_rubble <- chem_reframe_clean_RUBBLECCA %>%
  ggplot(aes(x = NEP_dailymean, y = log(deltaDOC_dailymean))) +
  geom_point(color = "tan") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEP" ~ (mmol ~ C ~ m^2 ~ h^-1))), y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1, label.y = 4.25, size = 5) + 
  stat_cor(label.x = 1, label.y = 4, size = 5)
NEP_deltaDOC_rubble

NEP_deltaDOC_rubble_model <- lmer(log(deltaDOC_dailymean) ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean_RUBBLECCA)
check_model(NEP_deltaDOC_rubble_model)
summary(NEP_deltaDOC_rubble_model)

# NEC v delta pH # 
NEC_deltapH_rubble <- chem_reframe_clean_RUBBLECCA %>%
  ggplot(aes(x = NEC_dailymean, y = deltapH_dailymean)) +
  geom_point(color = "tan") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "∆ pH") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 0.5, label.y = 0.05, size = 5) + 
  stat_cor(label.x = 0.5, label.y = 0.045, size = 5)
NEC_deltapH_rubble

NEC_deltapH_rubble_model <- lmer(deltapH_dailymean ~ NEC_dailymean + (1|TANK_NUM), data = chem_reframe_clean_RUBBLECCA)
check_model(NEC_deltapH_rubble_model)
summary(NEC_deltapH_rubble_model)

# NEC v delta DOC # 
NEC_deltaDOC_rubble <- chem_reframe_clean_RUBBLECCA %>%
  ggplot(aes(x = NEC_dailymean, y = log(deltaDOC_dailymean))) +
  geom_point(color = "tan") + 
  #geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1))), y = expression(bold("∆ DOC"))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 0.5, label.y = 4, size = 5) + 
  stat_cor(label.x = 0.5, label.y = 3.75, size = 5)
NEC_deltaDOC_rubble

NEC_deltaDOC_rubble_model <- lmer(log(deltaDOC_dailymean) ~ NEC_dailymean + (1|TANK_NUM), data = chem_reframe_clean_RUBBLECCA)
check_model(NEC_deltaDOC_rubble_model)
summary(NEC_deltaDOC_rubble_model)

## plot for RUBBLE/CCA ##
NEC_NEP_pH_DOC_RUBBLECCA <- (NEC_deltapH_rubble + NEP_deltapH_rubble)/(NEC_deltaDOC_rubble + NEP_deltaDOC_rubble) + plot_annotation(tag_levels = "a")
NEC_NEP_pH_DOC_RUBBLECCA
#ggsave(plot = NEC_NEP_pH_DOC_RUBBLECCA, filename = here("Output", "NEC_NEP_pH_DOC_RUBBLECCApatch.png"), width = 15, height = 15)

### NEC AND NEP V DELTA PH AND DELTA DOC DAY ###

## NEP v pH DAY ## 
NEP_deltapH_DAY_plot <- chem_reframe_DAY_clean %>%
  ggplot(aes(x = NEP, y = deltapH)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 1, label.y = 0.35, size = 5) + 
  stat_cor(label.x = 1, label.y = 0.3, size = 5) 
NEP_deltapH_DAY_plot
#ggsave(plot = NEP_deltapH_DAY_plot, filename = here("Output", "NEP_Plots", "NEP_deltapH_DAY_plot.png"), width = 10, height = 10)

## NEP v DOC DAY ##
NEP_deltaDOC_DAY_plot <- chem_reframe_DAY_clean %>%
  ggplot(aes(x = NEP, y = log(deltaDOC))) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEP" ~ (mmol ~ C ~ m^2 ~ h^-1))), y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 2, label.y = 1.75, size = 5) + 
  stat_cor(label.x = 2, label.y = 1.25, size = 5) 
NEP_deltaDOC_DAY_plot
#ggsave(plot = NEP_deltaDOC_DAY_plot, filename = here("Output", "NEP_Plots", "NEP_deltaDOC_DAY_plot.png"), width = 10, height = 10)

## NEC v pH DAY ## 
NEC_deltapH_DAY_plot <- chem_reframe_DAY_clean %>%
  ggplot(aes(x = NEC, y = deltapH)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "∆ pH") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 0, label.y = 0.35, size = 5) + 
  stat_cor(label.x = 0, label.y = 0.275, size = 5) 
NEC_deltapH_DAY_plot
#ggsave(plot = NEC_deltapH_DAY_plot, filename = here("Output", "TA_NECPlots", "NEC_deltapH_DAY_plot.png"), width = 10, height = 10)

## NEC v DOC DAY ##
NEC_deltaDOC_DAY_plot <- chem_reframe_DAY_clean %>%
  ggplot(aes(x = NEC, y = log(deltaDOC))) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEC" ~ mmol ~ CaCO[3] ~ m^2 ~ h^-1)), y = "∆ DOC") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = -1, label.y = 0.75, size = 5) + 
  stat_cor(label.x = -1, label.y = 0.4, size = 5) 
NEC_deltaDOC_DAY_plot
#ggsave(plot = NEC_deltaDOC_DAY_plot, filename = here("Output", "TA_NECPlots", "NEC_deltaDOC_DAY_plot.png"), width = 10, height = 10)

## NEC and NEP v delta pH and delta DOC daytime patchwork plot ##
NEC_NEP_pH_DOC_DAY <- (NEC_deltapH_DAY_plot + NEP_deltapH_DAY_plot)/(NEC_deltaDOC_DAY_plot + NEP_deltaDOC_DAY_plot) + plot_annotation(tag_levels = "a")
NEC_NEP_pH_DOC_DAY
#ggsave(plot = NEC_NEP_pH_DOC_DAY, filename = here("Output", "NEC_NEP_deltapH_DOC_DAY_patchwork.png"), width = 15, height = 15)

### NEC AND NEP V DELTA PH AND DELTA DOC NIGHT ###

## NEP v pH NIGHT ## 
NEP_deltapH_NIGHT_plot <- chem_reframe_NIGHT_clean %>%
  ggplot(aes(x = NEP, y = deltapH)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = -0.5, label.y = -0.15, size = 5) + 
  stat_cor(label.x = -0.5, label.y = -0.2, size = 5) 
NEP_deltapH_NIGHT_plot
#ggsave(plot = NEP_deltapH_NIGHT_plot, filename = here("Output", "NEP_Plots", "NEP_deltapH_NIGHT_plot.png"), width = 10, height = 10)

## NEP v DOC NIGHT ##
NEP_deltaDOC_NIGHT_plot <- chem_reframe_NIGHT_clean %>%
  ggplot(aes(x = NEP, y = log(deltaDOC))) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEP" ~ (mmol ~ C ~ m^2 ~ h^-1))), y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 0, label.y = 2, size = 5) + 
  stat_cor(label.x = 0, label.y = 1.75, size = 5) 
NEP_deltaDOC_NIGHT_plot
#ggsave(plot = NEP_deltaDOC_NIGHT_plot, filename = here("Output", "NEP_Plots", "NEP_deltaDOC_NIGHT_plot.png"), width = 10, height = 10)

## NEC v pH NIGHT ## 
NEC_deltapH_NIGHT_plot <- chem_reframe_NIGHT_clean %>%
  ggplot(aes(x = NEC, y = deltapH)) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "", y = "∆ pH") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 0.1, label.y = -0.15, size = 5) + 
  stat_cor(label.x = 0.1, label.y = -0.175, size = 5) 
NEC_deltapH_NIGHT_plot
#ggsave(plot = NEC_deltapH_NIGHT_plot, filename = here("Output", "TA_NECPlots", "NEC_deltapH_NIGHT_plot.png"), width = 10, height = 10)

## NEC v DOC NIGHT ##
NEC_deltaDOC_NIGHT_plot <- chem_reframe_NIGHT_clean %>%
  ggplot(aes(x = NEC, y = log(deltaDOC))) +
  geom_point(aes(color = TREATMENT)) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) +
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = expression(bold("NEC" ~ mmol ~ CaCO[3] ~ m^2 ~ h^-1)), y = "∆ DOC") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none") +
  stat_regline_equation(label.x = 0.1, label.y = 2, size = 5) + 
  stat_cor(label.x = 0.1, label.y = 1.75, size = 5) 
NEC_deltaDOC_NIGHT_plot
#ggsave(plot = NEC_deltaDOC_NIGHT_plot, filename = here("Output", "TA_NECPlots", "NEC_deltaDOC_NIGHT_plot.png"), width = 10, height = 10)

## NEC and NEP v delta pH and delta DOC nighttime patchwork plot ##
NEC_NEP_pH_DOC_NIGHT <- (NEC_deltapH_NIGHT_plot + NEP_deltapH_NIGHT_plot)/(NEC_deltaDOC_NIGHT_plot + NEP_deltaDOC_NIGHT_plot) + plot_annotation(tag_levels = "a")
NEC_NEP_pH_DOC_NIGHT
#ggsave(plot = NEC_NEP_pH_DOC_NIGHT, filename = here("Output", "NEC_NEP_deltapH_DOC_NIGHT_patchwork.png"), width = 15, height = 15)


