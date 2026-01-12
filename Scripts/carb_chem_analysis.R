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
library(pwr)
#install.packages("pwr")
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


## write full chemistry data file ##
# temp, flow, salinity, pH, light, DIC, pH inflow, TA inflow, DIC inflow, 
# delta pH, delta TA, delta DIC, NEC, NEP, DOC 

raw_chem_data <- Data %>%
  filter(TIME %in% c("12:00:00", "21:00:00")) %>%
  group_by(TREATMENT, TANK_NUM) %>% 
  summarize(grand_mean_pH = mean(pH, na.rm = TRUE),
            grand_mean_DOC = mean(NPOC_uM, na.rm = TRUE),
            grand_range_pH = max(pH, na.rm = TRUE) - min(pH, na.rm = TRUE), 
            grand_range_DOC = max(NPOC_uM, na.rm = TRUE) - min(NPOC_uM, na.rm = TRUE),
            grand_max_pH = max(pH, na.rm = TRUE), 
            grand_max_DOC = max(NPOC_uM, na.rm = TRUE),
            grand_min_pH = min(pH, na.rm = TRUE), 
            grand_min_DOC = min(NPOC_uM, na.rm = TRUE))%>%
  drop_na()

write_csv(raw_chem_data, here("Data", "Chemistry", "Raw_Chem_Data.csv"))

write_csv(Data, here("Data", "Chemistry", "Full_Carb_Chem_Data.csv"))
full_carb_chem_data <- read_csv(here("Data", "Chemistry", "Full_Carb_Chem_Data.csv"))
full_carb_chem_data$INFLOW_TABLE <- as.factor(full_carb_chem_data$INFLOW_TABLE)

tank_residence_time <- Data %>%
  group_by(TANK_NUM) %>%
  summarize(avg_res_time = mean(residence_time, na.rm = TRUE),
            res_error = sd(residence_time, na.rm = TRUE)/sqrt(n()),
            avg_flow_rate = mean(flowrate, na.rm = TRUE),
            flow_error = sd(flowrate, na.rm = TRUE)/sqrt(n())) %>%
  drop_na()
tank_residence_time
write_csv(tank_residence_time, here("Data", "Chemistry", "Tank_Residence_Time.csv"))

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
         NEC = ifelse(NEC < -1, NA, NEC),
         NEP = ifelse(NEP > 2, NA, NEP), 
         NEP = ifelse(TREATMENT == "Control" & NEP > 1, NA, NEP),
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

        
#write_csv(chem_summary_data, here("Data", "Chemistry", "chem_summary_data.csv"))

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

########## create geom rect ############
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

###################################
### NEC DATA ### 
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
  labs(x = "", y = expression(bold("NCC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Macroalgae-Enriched", "Control" = "Control",
                            "Coral_Dom" = "Coral-Enriched", "Rubble_Dom" = "CCA-Enriched")) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = NEC_data_day, aes(x = TREATMENT, y = NEC_day_mean), size = 3, alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 5) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")
NEC_day_mean_plot
ggsave(plot = NEC_day_mean_plot, filename = here("Output", "TA_NECPlots", "NEC_day_mean_plot.png"), width = 9, height = 9)

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
  labs(x = "", y = "") +
  scale_x_discrete(labels=c("Algae_Dom" = "Macroalgae-Enriched", "Control" = "Control",
                            "Coral_Dom" = "Coral-Enriched", "Rubble_Dom" = "CCA-Enriched")) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = NEC_data_night, aes(x = TREATMENT, y = NEC_night_mean), size = 3, alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 5) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")
NEC_night_mean_plot
ggsave(plot = NEC_night_mean_plot, filename = here("Output", "TA_NECPlots", "NEC_night_mean_plot.png"), width = 9, height = 9)

NEC_data_night$TREATMENT <- factor(NEC_data_night$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEC_nighttime_model <- lmer(NEC_night_mean ~ TREATMENT + (1|TANK_NUM), data=NEC_data_night)
check_model(NEC_nighttime_model)
summary(NEC_nighttime_model)
anova(NEC_nighttime_model) 

emmeans(NEC_nighttime_model, pairwise ~ "TREATMENT", adjust = "Tukey")
NEC_nighttime_model_noRandom <- lm(NEC_night_mean ~ TREATMENT, data = NEC_data_night)
HSD.test(NEC_nighttime_model_noRandom, "TREATMENT", console=TRUE)

##### NEP DATA ANALYSIS ##### 
#############NEP #######################

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
  ggtitle("Daytime Mean") +
  geom_hline(yintercept = 0, lty =2) +
  ylim(c(-3,5)) +
  labs(x = "", y = expression(bold("NCP" ~ (mmol ~ C ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Macroalgae-Enriched", "Control" = "Control",
                            "Coral_Dom" = "Coral-Enriched", "Rubble_Dom" = "CCA-Enriched")) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = NEP_data_day, aes(x = TREATMENT, y = NEP_day_mean), size = 3, alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 5) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, hjust = 0.5), 
        legend.position = "none")
NEP_day_mean_plot
ggsave(plot = NEP_day_mean_plot, filename = here("Output", "NEP_Plots", "NEP_day_mean_plot.png"), width = 9, height = 9)

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
  ggtitle("Nighttime Mean") +
  geom_hline(yintercept = 0, lty = 2) +
  ylim(c(-2,2)) +
  labs(x = "", y = "") +
  scale_x_discrete(labels=c("Algae_Dom" = "Macroalgae-Enriched", "Control" = "Control",
                            "Coral_Dom" = "Coral-Enriched", "Rubble_Dom" = "CCA-Enriched")) +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan")) +
  geom_point(data = NEP_data_night, aes(x = TREATMENT, y = NEP_night_mean), size = 3, alpha = 0.25) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.25) +
  stat_summary(fun.y = mean, geom = "point", size = 5) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none")
NEP_night_mean_plot
ggsave(plot = NEP_night_mean_plot, filename = here("Output", "NEP_Plots", "NEP_night_mean_plot.png"), width = 9, height = 12)

NEP_data_night$TREATMENT <- factor(NEP_data_night$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEP_night_mean_model <- lmer(NEP_night_mean ~ TREATMENT + (1|TANK_NUM), data=NEP_data_night)
check_model(NEP_night_mean_model)
summary(NEP_night_mean_model)
anova(NEP_night_mean_model) # non sig effect of community on nighttime nep. p = 0.2065


## create patchwork plot
NEP_v_NEC_patchwork <- (NEP_day_mean_plot + NEP_night_mean_plot)/(NEC_day_mean_plot + NEC_night_mean_plot) + plot_annotation(tag_levels = "a")
NEP_v_NEC_patchwork
ggsave(plot = NEP_v_NEC_patchwork, filename = here("Output", "NEP_v_NEC_patchwork.png"), width = 12, height = 10)

############# NEP DATA ############

#### Look at NEP vs pH and DOC
Clean_Chem_all<- chem_reframe_NIGHT_clean %>% 
  bind_rows(chem_reframe_DAY_clean) %>%
  filter(NEP < 4)

# ANCOVA with pH~NEP*Treament and accounting for tankID

Clean_Chem_all2 <- Clean_Chem_all %>%
  mutate(DATE = ymd(paste(DATE)))

clean_chem_all_sum <- Clean_Chem_all %>%
  group_by(TREATMENT) %>%
  summarise(Count = n())
clean_chem_all_sum
pwr.anova.test(k = 4, n = 28, f = 0.25, sig.level = 0.05, power = NULL)

mod_DOC<-lm(DOC~NEP*TREATMENT + (1|DATE), data = Clean_Chem_all)
anova(mod_DOC)
summary(mod_DOC)

mod_DOC_NEC<-lm(DOC~NEC*TREATMENT + (1|DATE), data = Clean_Chem_all)
anova(mod_DOC_NEC)
summary(mod_DOC_NEC)

mod_NEP_NEC<-lm(NEC ~ pH*TREATMENT + (1|DATE), data = Clean_Chem_all)
anova(mod_NEP_NEC)
summary(mod_NEP_NEC)

Clean_Chem_all$TREATMENT <- factor(Clean_Chem_all$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

chem_plot_all_pH<-Clean_Chem_all %>%
  ggplot(aes(x = NEP, y = pH))+
  geom_point(aes(color = TREATMENT), size = 3, alpha = 0.5)+
  geom_smooth(method = "lm", color = "black")+
  labs(color = "",
       y = expression("pH"[T]),
       x = expression("NCP (mmol C m"^2~"hr"^-1~")"))+
  theme_bw() + 
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none") +
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched",
                                "CCA-Enriched"), values = c("blue", "darkgreen", "coral", "tan"), guide = FALSE)
chem_plot_all_pH


chem_plot_all_DOC <- Clean_Chem_all %>%
  ggplot(aes(x = NEP, y = DOC))+
  geom_point(aes(color = TREATMENT), size = 3, alpha = 0.5)+
  geom_smooth(method = "lm", color = "black")+
  labs(color = "",
       x = expression("NCP (mmol C m"^2~"hr"^-1~")"),
       y = expression("DOC ("~mu~"mol L"^-1~")")) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 14)) +
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched",
                                "CCA-Enriched"), values = c("blue", "darkgreen", "coral", "tan"))

chem_plot_all_DOC


chem_plot_all_pH_NEC<-Clean_Chem_all %>%
  ggplot(aes(x = pH, y = NEC))+
  geom_point(aes(color = TREATMENT), size = 3, alpha = 0.5)+
  geom_smooth(method = "lm", color = "black")+
  labs(color = "",
       y = expression("NCC (mmol C m"^2~" hr"^-1~")"),
       x = expression("pH"[T]))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none") +
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched",
                                "CCA-Enriched"), values = c("blue", "darkgreen", "coral", "tan"), guide = FALSE)
chem_plot_all_pH_NEC

chem_plot_all_DOC_NEC<-Clean_Chem_all %>%
  ggplot(aes(x = NEC, y = DOC))+
  geom_point(aes(color = TREATMENT), size = 3, alpha = 0.5)+
  labs(color = "",
       y = expression("DOC ("~mu~"mol L"^-1~")"),
       x = expression("NCC (mmol C m"^2~" hr"^-1~")"))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 14)) +
  scale_color_manual(labels = c("Control", "Macroalgae-Enriched", "Coral-Enriched",
                                "CCA-Enriched"), values = c("blue", "darkgreen", "coral", "tan"))
chem_plot_all_DOC_NEC

ncp_pH_DOC_patch <- (chem_plot_all_pH+chem_plot_all_pH_NEC)/(chem_plot_all_DOC+chem_plot_all_DOC_NEC)  +
  plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = 'bottom') 
ncp_pH_DOC_patch

ggsave(plot = ncp_pH_DOC_patch, filename = here("Output", "ncp_pH_DOC_patch.png"), width = 12, height = 12)