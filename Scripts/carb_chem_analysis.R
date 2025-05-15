library(tidyverse)
library(seacarb)
library(broom)
library(here)
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
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  right_join(.,pHData) %>% # join with the pH sample data
  mutate(mVTris = TEMPINLAB*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  drop_na(TEMPINSITU)%>%
  drop_na(mV) %>%
  mutate(pH = pH(Ex=mV,Etris=mVTris,S=SALINITY,T=TEMPINLAB)) %>%  # calculate pH of the samples using the pH seacarb function
#Now calculate pH
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
  mutate(DIC_mmol_kg = (DIC*1e6)) 

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
         DIC_inflow = DIC_mmol_kg,
         DOC_inflow = NPOC_uM) %>%
  mutate(INFLOW_TABLE = ifelse(TANK_NUM == "Inflow1",1,2)) %>% # give them inflow numbers to pair easily with the TankID 
  ungroup()%>%
  select(DATE,TIME, INFLOW_TABLE, pH_inflow, TA_inflow, DIC_inflow, DOC_inflow) # drop the Tank ID column to be able to join the data correctly by inflow #

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
         deltaTA = TA_inflow - TA, # calculate the difference between in and outflow
         deltaDIC = DIC_inflow - DIC_mmol_kg, # calculate delta DIC to calculate NEP 
         NEC = (deltaTA/2)*(1.025)*(10)*(1/residence_time)*(1/SurfaceArea), ### for a real rate should probably normalize the delta TA to the delta control just like in respo
         NEP = ((deltaDIC)*(1.025)*(10)*(1/residence_time)*(1/SurfaceArea)) - NEC
         )

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


## create one summary data frame with means and ranges for TA, pH, and DOC ##
# start by isolating 12:00 and 21:00 time periods #
## filtering for 12:00 and 21:00 sampling times and reframing to add daily mean and daily range between 12 and 9 

chem_reframe_date_removed <- Data %>% # chem reframe with 20240606 removed from analyses 
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  filter(!DATE == "20240606") %>%
  group_by(TREATMENT, DATE, TANK_NUM) %>%
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

chem_reframe_date_remove_plot <- chem_reframe_date_removed %>%
  select(TREATMENT, TANK_NUM, TA_range:NEP_dailymean) %>%
  pivot_longer(cols = TA_range:NEP_dailymean) %>%
  ggplot(aes(x = TREATMENT, y = value)) +
  geom_point() + 
  facet_wrap(~name, scales = "free")
chem_reframe_date_remove_plot

chem_reframe <- Data %>% 
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
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

chem_reframe_plots <- chem_reframe %>%
  select(TREATMENT, TANK_NUM, TA_range:NEP_dailymean) %>%
  pivot_longer(cols = TA_range:NEP_dailymean) %>%
  ggplot(aes(x = TREATMENT, y = value)) +
  geom_point() + 
  facet_wrap(~name, scales = "free")
chem_reframe_plots

chem_reframe_clean <- chem_reframe %>%
  mutate(deltaTA_range = ifelse(deltaTA_range < -250, NA, deltaTA_range),
         deltaTA_dailymean = ifelse(deltaTA_dailymean < -250, NA, deltaTA_dailymean),
         deltapH_dailymean = ifelse(deltapH_dailymean > 0.4, NA, deltapH_dailymean), 
         deltapH_range = ifelse(deltapH_range < -0.15, NA, deltapH_range), 
         NEC_dailymean = ifelse(NEC_dailymean < -4, NA, NEC_dailymean), 
         NEC_dailymean = ifelse(TREATMENT ==  "Coral_Dom"&NEC_dailymean < -1, NA, NEC_dailymean),
         NEC_range = ifelse(NEC_range < -4, NA, NEC_range),
         NEP_range = ifelse(NEP_range < -6, NA, NEP_range), 
         NEP_dailymean = ifelse(NEP_dailymean < -2.5, NA, NEP_dailymean),
         TA_dailymean = ifelse(TA_dailymean > 2500, NA, TA_dailymean), 
         TA_dailymean = ifelse(TREATMENT == "Coral_Dom"&TA_dailymean < 2200, NA, TA_dailymean),
         TA_range = ifelse(TA_range > 400, NA, TA_range))

chem_reframe_date_removed_clean <- chem_reframe_date_removed %>%
  mutate(DOC_dailymean = ifelse(TREATMENT == "Rubble_Dom"&DOC_dailymean > 300, NA, DOC_dailymean),
         NEC_range = ifelse(NEC_range < -0.4, NA, NEC_range),
         NEP_range = ifelse(NEP_range < -3, NA, NEP_range), 
         TA_dailymean = ifelse(TREATMENT == "Coral_Dom"&TA_dailymean < 2200, NA, TA_dailymean),
         TA_range = ifelse(TA_range < -150, NA, TA_range))


chem_reframe_clean %>%
  select(TREATMENT, TANK_NUM, TA_range:NEP_dailymean) %>%
  pivot_longer(cols = TA_range:NEP_dailymean) %>%
  ggplot(aes(x = TREATMENT, y = value)) +
  geom_point() + 
  facet_wrap(~name, scales = "free")


chem_reframe_date_removed_clean %>%
  select(TREATMENT, TANK_NUM, TA_range:NEP_dailymean) %>%
  pivot_longer(cols = TA_range:NEP_dailymean) %>%
  ggplot(aes(x = TREATMENT, y = value)) +
  geom_point() + 
  facet_wrap(~name, scales = "free")


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

ggplot(chem_reframe_clean) + 
  geom_point(aes(x = NEP_dailymean, y = NEC_dailymean, color = TREATMENT)) + 
  facet_wrap(~TREATMENT, scales = "free")

ggplot(chem_reframe_date_removed_clean) + 
  geom_point(aes(x = NEP_dailymean, y = NEC_dailymean, color = TREATMENT)) + 
  facet_wrap(~TREATMENT, scales = "free")


ggplot(chem_reframe_clean) + 
  geom_point(aes(x = NEP_dailymean, y = pH_dailymean, color = TREATMENT)) + 
  facet_wrap(~TREATMENT, scales = "free")

NEP_pH_model <- lmer(pH_dailymean ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean)
check_model(NEP_pH_model)
summary(NEP_pH_model)


chem_reframe_date_removed_clean$TREATMENT <- factor(chem_reframe_date_removed_clean$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"), 
                                                    labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated"))

NEP_pH_plot <- ggplot(chem_reframe_date_removed_clean) + 
  labs(x = "Daily Mean NEP", y = "Daily Mean pH") +
  geom_point(aes(x = NEP_dailymean, y = pH_dailymean, color = TREATMENT)) + 
  geom_smooth(aes(x = NEP_dailymean, y = pH_dailymean)) + 
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"))
NEP_pH_plot

#ggsave(plot = NEP_pH_plot, filename = here("Output", "NEP_plots", "NEP_pH_plot.png"), width = 9, height = 6)

NEP_pH_model2 <- lmer(pH_dailymean ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_date_removed_clean)
check_model(NEP_pH_model2)
summary(NEP_pH_model2)


ggplot(chem_reframe_clean) + 
  geom_point(aes(x = NEC_dailymean, y = TA_dailymean, color = TREATMENT)) + 
  facet_wrap(~TREATMENT, scales = "free")

NEC_TA_plot <- ggplot(chem_reframe_date_removed_clean) + 
  labs(x = "Daily Mean NEC", y = "Daily Mean TA") +
  geom_point(aes(x = NEC_dailymean, y = TA_dailymean, color = TREATMENT)) + 
  geom_smooth(aes(x = NEC_dailymean, y = TA_dailymean)) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"))
NEC_TA_plot
#ggsave(plot = NEC_TA_plot, filename = here("Output", "TA_NECPlots", "NEC_TA_plot.png"), width = 9, height = 6)


NEC_TA_model <- lmer(TA_dailymean ~ NEC_dailymean * TREATMENT + (1|TANK_NUM), data = chem_reframe_clean)
check_model(NEC_TA_model)
summary(NEC_TA_model)

NEC_TA_model2 <- lmer(TA_dailymean ~ NEC_dailymean + (1|TANK_NUM), data = chem_reframe_date_removed_clean)
check_model(NEC_TA_model2)
summary(NEC_TA_model2)



NEP_DOC_plot <- ggplot(chem_reframe_date_removed_clean) + 
  labs(x = "Daily Mean NEP", y = "Daily Mean DOC") +
  geom_point(aes(x = NEP_dailymean, y = DOC_dailymean, color = TREATMENT)) + 
  geom_smooth(aes(x = NEP_dailymean, y = DOC_dailymean)) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) + 
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan")) 
NEP_DOC_plot
#ggsave(plot = NEP_DOC_plot, filename = here("Output", "NEP_Plots", "NEP_DOC_plot.png"), width = 9, height = 6)

NEP_DOC_model <- lmer(DOC_dailymean ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_clean)
check_model(NEP_DOC_model)
summary(NEP_DOC_model)

NEP_DOC_model2 <- lmer(DOC_dailymean ~ NEP_dailymean + (1|TANK_NUM), data = chem_reframe_date_removed_clean)
check_model(NEP_DOC_model2)
summary(NEP_DOC_model2)


ggplot(chem_reframe_date_removed_clean) + 
  geom_point(aes(x = NEP_dailymean, y = DOC_dailymean, color = TREATMENT)) + 
  facet_wrap(~TREATMENT, scales = "free")


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

NEC_plot <- chem_summary_data %>% # chem reframe data is only 12 and 9 pm sampling, CLEANED 
  ggplot(aes(x = TREATMENT, y = NEC_mean, color = TREATMENT)) +
  theme_classic() +
  labs(x="Dominant Benthic Community",
       y = expression(bold("Daily Mean Net Ecosystem Calcification (NEC)" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels = c("Control", "Algae-Dominated", "Coral-Dominated", "Rubble/CCA-Dominated")) +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) + 
  geom_jitter(data = chem_reframe_clean, aes(x = TREATMENT, y = NEC_dailymean), alpha = 0.7) +
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1, color = "black") +
  scale_color_manual(values = c("Control" = "blue", "Algae_Dom" = "darkgreen", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEC_plot
#ggsave(plot = NEC_plot, filename = here("Output", "TA_NECPlots", "NEC_mean_plot.png"), width = 9, height = 10)

## NEC daily mean stats ##
NEC_mean_model <- lmer(NEC_dailymean ~ TREATMENT + (1|TANK_NUM), data=chem_reframe_clean)
check_model(NEC_mean_model)
summary(NEC_mean_model) # coral and rubble/cca community significantly different from algae 
anova(NEC_mean_model) # significant effect of community type on daily mean NEC 
emmeans(NEC_mean_model, pairwise ~ "TREATMENT", adjust = "Tukey")

NEC_mean_model_noRandom <- lm(NEC_dailymean ~ TREATMENT, data = chem_reframe_clean)
HSD.test(NEC_mean_model_noRandom, "TREATMENT", console=TRUE)

## NEC range ## 
## plot daily NEC range ## 
NEC_range_plot <- chem_summary_data %>%
  ggplot(aes(x = TREATMENT, y = NEC_rangemean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daily Mean NEC Range" ~ (mmol ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble/CCA-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = chem_reframe_clean, aes(x = TREATMENT, y = NEC_range), alpha = 0.7) +
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width = 0.1, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEC_range_plot 
#ggsave(plot = NEC_range_plot, filename = here("Output", "NEC_range_plot.png"), width = 9, height = 6)

NEC_range_model <- lmer(NEC_range ~ TREATMENT + (1|TANK_NUM), data=chem_reframe_clean)
check_model(NEC_range_model)
summary(NEC_range_model)
anova(NEC_range_model) 




## NEC DAY vs NEC NIGHT ##
# NEC day

NEC_data_day <- Data %>%
  filter(TIME %in% "12:00:00") %>%
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  select(DATE, TIME, DATETIME, TANK_NUM, TREATMENT, NEC) %>%
  filter(!(NEC < -2 & TREATMENT == "Coral_Dom")) %>%
  reframe(NEC_day_mean = mean(NEC, na.rm = TRUE)) %>%
  drop_na()

NEC_day_plotdata <- NEC_data_day %>%
  group_by(TREATMENT) %>%
  summarize(NEC_mean = mean(NEC_day_mean, na.rm = TRUE),
            NEC_se = sd(NEC_day_mean, na.rm = TRUE)/sqrt(n()))

NEC_day_plotdata$TREATMENT <- factor(NEC_day_plotdata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEC_day_mean_plot <- NEC_day_plotdata %>%
  ggplot(aes(x = TREATMENT, y = NEC_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daytime Mean NEC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = NEC_data_day, aes(x = TREATMENT, y = NEC_day_mean), alpha = 0.7) +
  geom_errorbar(aes(ymin = NEC_mean - NEC_se,
                    ymax = NEC_mean + NEC_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEC_day_mean_plot
#ggsave(plot = NEC_day_mean_plot, filename = here("Output", "TA_NECPlots", "NEC_day_mean_plot.png"), width = 9, height = 9)

NEC_daytime_model <- lmer(NEC_day_mean ~ TREATMENT + (1|TANK_NUM), data=NEC_data_day)
check_model(NEC_daytime_model)
summary(NEC_daytime_model)
anova(NEC_daytime_model) # sig effect of community on daytime NEC rates

emmeans(NEC_daytime_model, pairwise ~ "TREATMENT", adjust = "Tukey")

NEC_daytime_model_noRandom <- lm(NEC_day_mean ~ TREATMENT, data = NEC_data_day)
HSD.test(NEC_daytime_model_noRandom, "TREATMENT", console=TRUE)


# NEC night
NEC_data_night <- Data %>%
  filter(TIME %in% "21:00:00") %>%
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  select(DATE, TIME, DATETIME, TANK_NUM, TREATMENT, NEC) %>%
  filter(!(NEC < -2 & TREATMENT == "Coral_Dom")) %>%
  reframe(NEC_night_mean = mean(NEC, na.rm = TRUE)) %>%
  drop_na()

NEC_night_plotdata <- NEC_data_night %>%
  group_by(TREATMENT) %>%
  summarize(NEC_mean_n = mean(NEC_night_mean, na.rm = TRUE),
            NEC_se_n = sd(NEC_night_mean, na.rm = TRUE)/sqrt(n()))

NEC_night_plotdata$TREATMENT <- factor(NEC_night_plotdata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEC_night_mean_plot <- NEC_night_plotdata %>%
  ggplot(aes(x = TREATMENT, y = NEC_mean_n, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Nighttime Mean NEC" ~ (mmol ~ CaCO[3] ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = NEC_data_night, aes(x = TREATMENT, y = NEC_night_mean), alpha = 0.7) +
  geom_errorbar(aes(ymin = NEC_mean_n - NEC_se_n,
                    ymax = NEC_mean_n + NEC_se_n), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEC_night_mean_plot
#ggsave(plot = NEC_night_mean_plot, filename = here("Output", "TA_NECPlots", "NEC_night_mean_plot.png"), width = 9, height = 9)

NEC_data_night$TREATMENT <- factor(NEC_data_night$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEC_nighttime_model <- lmer(NEC_night_mean ~ TREATMENT + (1|TANK_NUM), data=NEC_data_night)
check_model(NEC_nighttime_model)
summary(NEC_nighttime_model)
anova(NEC_nighttime_model)


##### NEP DATA ANALYSIS ##### 

NEP_data <- Data %>%
  select(DATETIME, DATE, TIME, NEP, TREATMENT, TANK_NUM) %>%
  group_by(DATETIME, DATE, TIME, TREATMENT, TANK_NUM) %>%
  drop_na()

NEP_data2 <- NEP_data %>% 
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(NEP_range = NEP[TIME == hms("12:00:00")] - NEP[TIME == hms("21:00:00")],
          NEP_dailymean = mean(NEP, na.rm = TRUE)) %>%
  filter(NEP_range > -5) %>%
  filter(!(TREATMENT == "Coral_Dom" & NEP_range < -2)) %>%
  filter(!NEP_range < -2)

#create plot data
NEP_plotdata <- NEP_data2 %>%
  group_by(TREATMENT) %>%
  summarize(NEP_rangemean = mean(NEP_range, na.rm = TRUE),
            NEP_rangese = sd(NEP_range, na.rm = TRUE)/sqrt(n()),
            NEP_mean = mean(NEP_dailymean, na.rm = TRUE),
            NEP_se = sd(NEP_dailymean, na.rm = TRUE)/sqrt(n())) 
NEP_plotdata

NEP_plotdata2 <- NEP_data2 %>%
  group_by(TREATMENT, TANK_NUM) %>%
  summarize(NEP_rangemean = mean(NEP_range, na.rm = TRUE), 
            NEP_rangese = sd(NEP_range, na.rm = TRUE)/sqrt(n()),
            NEP_mean = mean(NEP_dailymean, na.rm = TRUE), 
            NEP_se = sd(NEP_dailymean, na.rm = TRUE)/sqrt(n()))
NEP_plotdata2

# NEP daily range plot
NEP_plotdata$TREATMENT <- factor(NEP_plotdata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

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
anova(NEP_range_model)

emmeans(NEP_range_model, pairwise ~ "TREATMENT", adjust = "Tukey")

NEP_range_model_noRandom <- lm(NEP_range ~ TREATMENT, data = NEP_data2)
HSD.test(NEP_range_model_noRandom, "TREATMENT", console=TRUE)

# NEP mean #
NEP_mean_plot <- NEP_plotdata %>%
  ggplot(aes(x = TREATMENT, y = NEP_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daily Mean NEP" ~ (mmol ~ C ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble/CCA-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = NEP_data2, aes(x = TREATMENT, y = NEP_dailymean), alpha = 0.7) +
  geom_errorbar(aes(ymin = NEP_mean - NEP_se,
                    ymax = NEP_mean + NEP_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEP_mean_plot

#ggsave(plot = NEP_mean_plot, filename = here("Output", "NEP_Plots", "NEP_mean_plot.png"), width = 9, height = 9)

NEP_mean_model <- lmer(NEP_dailymean ~ TREATMENT + (1|TANK_NUM), data=NEP_data2)
check_model(NEP_mean_model)
summary(NEP_mean_model)
anova(NEP_mean_model)

## NEP DAY TIME ## 
NEP_data_day <- NEP_data %>% 
  filter(TIME %in% "12:00:00") %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(NEP_day_mean = mean(NEP, na.rm = TRUE)) %>%
  filter(!NEP_day_mean < -5) %>%
  filter(!(TREATMENT == "Coral_Dom" & NEP_day_mean < -1))

NEP_day_plotdata <- NEP_data_day %>%
  group_by(TREATMENT) %>%
  summarize(NEP_mean = mean(NEP_day_mean, na.rm = TRUE),
            NEP_se = sd(NEP_day_mean, na.rm = TRUE)/sqrt(n()))

NEP_day_plotdata$TREATMENT <- factor(NEP_day_plotdata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEP_day_mean_plot <- NEP_day_plotdata %>%
  ggplot(aes(x = TREATMENT, y = NEP_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Daytime Mean NEP" ~ (mmol ~ C ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble/CCA-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = NEP_data_day, aes(x = TREATMENT, y = NEP_day_mean), alpha = 0.7) +
  geom_errorbar(aes(ymin = NEP_mean - NEP_se,
                    ymax = NEP_mean + NEP_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEP_day_mean_plot

#ggsave(plot = NEP_day_mean_plot, filename = here("Output", "NEP_Plots", "NEP_day_mean_plot.png"), width = 9, height = 9)

NEP_data_day$TREATMENT <- factor(NEP_data_day$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEP_day_mean_model <- lmer(NEP_day_mean ~ TREATMENT + (1|TANK_NUM), data=NEP_data_day)
check_model(NEP_day_mean_model)
summary(NEP_day_mean_model)
anova(NEP_day_mean_model)

emmeans(NEP_day_mean_model, pairwise ~ "TREATMENT", adjust = "Tukey")

NEP_day_mean_model_noRandom <- lm(NEP_day_mean ~ TREATMENT, data = NEP_data_day)
HSD.test(NEP_day_mean_model_noRandom, "TREATMENT", console=TRUE)




## NEP NIGHT TIME ## 
NEP_data_night <- NEP_data %>%
  filter(TIME %in% "21:00:00") %>% 
  group_by(TREATMENT, DATE, TANK_NUM) %>%
  reframe(NEP_night_mean = mean(NEP, na.rm = TRUE))

NEP_night_plotdata <- NEP_data_night %>%
  group_by(TREATMENT) %>%
  summarize(NEP_mean = mean(NEP_night_mean, na.rm = TRUE),
            NEP_se = sd(NEP_night_mean, na.rm = TRUE)/sqrt(n()))

NEP_night_plotdata$TREATMENT <- factor(NEP_night_plotdata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEP_night_mean_plot <- NEP_night_plotdata %>%
  ggplot(aes(x = TREATMENT, y = NEP_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = expression(bold("Nighttime Mean NEP" ~ (mmol ~ C ~ m^2 ~ h^-1)))) +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble/CCA-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = NEP_data_night, aes(x = TREATMENT, y = NEP_night_mean), alpha = 0.7) +
  geom_errorbar(aes(ymin = NEP_mean - NEP_se,
                    ymax = NEP_mean + NEP_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
NEP_night_mean_plot

#ggsave(plot = NEP_night_mean_plot, filename = here("Output", "NEP_Plots", "NEP_night_mean_plot.png"), width = 9, height = 9)

NEP_data_night$TREATMENT <- factor(NEP_data_night$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

NEP_night_mean_model <- lmer(NEP_night_mean ~ TREATMENT + (1|TANK_NUM), data=NEP_data_night)
check_model(NEP_night_mean_model)
summary(NEP_night_mean_model)
anova(NEP_night_mean_model)

######################################
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

write_csv(pH_plotdata_tank_full, here("Data", "Chemistry", "Cleaned_pH_Data_per_Tank.csv"))


## plot pH range from 12:00 and 21:00 sampling throughout experiment ## 
pH_plotdata$TREATMENT <- factor(pH_plotdata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

pH_range_plot <- pH_plotdata %>%
  ggplot(aes(x = TREATMENT, y = pH_rangemean, color = TREATMENT)) +
  labs(x = "Treatment", y = "Daily Mean pH Range") +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = pH_clean, aes(x = TREATMENT, y = pH_range), alpha = 0.7) +
  geom_errorbar(aes(ymin = pH_rangemean - pH_rangese,
                    ymax = pH_rangemean + pH_rangese), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
pH_range_plot
ggsave(plot = pH_range_plot, filename = here("Output", "pH_range_plot.png"), width = 9, height = 6)

## pH range stats ##
pH_range_model <- lmer(pH_range ~ TREATMENT +(1|TANK_NUM), data=pH_clean)
check_model(pH_range_model)
summary(pH_range_model) # significant effect of control treatment on daily mean pH but not what we were 
# looking for. expected to see differences between the communities... 
anova(pH_range_model)

# take outlier from rubble dom community out and replot/redo stats 
pH_clean_filtered <- pH_clean %>%
  filter(pH_range > -0.2)

pH_plotdata2 <- pH_clean_filtered %>%
  group_by(TREATMENT) %>%
  summarize(pH_rangemean = mean(pH_range, na.rm = TRUE),
            pH_rangese = sd(pH_range, na.rm = TRUE)/sqrt(n()),
            pH_mean = mean(pH_dailymean, na.rm = TRUE),
            pH_se = sd(pH_dailymean, na.rm = TRUE)/sqrt(n()))
pH_plotdata2

pH_plotdata2$TREATMENT <- factor(pH_plotdata2$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

pH_range_plot2 <- pH_plotdata2 %>%
  ggplot(aes(x = TREATMENT, y = pH_rangemean, color = TREATMENT)) +
  labs(x = "Treatment", y = "Daily Mean pH Range") +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = pH_clean_filtered, aes(x = TREATMENT, y = pH_range), alpha = 0.7) +
  geom_errorbar(aes(ymin = pH_rangemean - pH_rangese,
                    ymax = pH_rangemean + pH_rangese), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") +
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
pH_range_plot2

# stats with filtered data set
pH_range_model2 <- lmer(pH_range ~ TREATMENT +(1|TANK_NUM), data=pH_clean_filtered)
check_model(pH_range_model2)
summary(pH_range_model2)
anova(pH_range_model2) # control still only the significant treatment on mean pH range 

## mean pH plot of 12:00 and 21:00 sampling throughout experiment ## 
pH_mean_plot <- pH_plotdata2 %>%
  ggplot(aes(x = TREATMENT, y = pH_mean, color = TREATMENT)) +
  labs(x = "Treatment", y = "Daily Mean pH") +
  scale_x_discrete(labels=c("Algae_Dom" = "Algae-Dominated", "Control" = "Control",
                            "Coral_Dom" = "Coral-Dominated", "Rubble_Dom" = "Rubble-Dominated")) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_jitter(data = pH_clean_filtered, aes(x = TREATMENT, y = pH_dailymean), alpha = 0.7) +
  geom_errorbar(aes(ymin = pH_mean - pH_se,
                    ymax = pH_mean + pH_se), color = "black", width = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2.5, color = "black") + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
pH_mean_plot
ggsave(plot = pH_mean_plot, filename = here("Output", "pH_mean_plot.png"), width = 9, height = 6)

## mean pH treatment stats ## 
mean_pH_model <- lmer(pH_dailymean~TREATMENT +(1|TANK_NUM), data=pH_clean_filtered)
check_model(mean_pH_model)
summary(mean_pH_model)
anova(mean_pH_model)
# no significant effect of community type on daily mean pH 


##### NEP and pH #####
Data1 <- Data %>%
  filter(!NEP < -4) %>%
  filter(!TA > 2700)

Data1$TREATMENT <- factor(Data1$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))
NEP_pH_plot <- Data1 %>%
  ggplot(aes(x = NEP, y = pH)) +
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "NEP", y = "pH") +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 0, label.y = 8.45, size = 5) + 
  stat_cor(label.x = 0, label.y = 8.4, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"))
NEP_pH_plot
#ggsave(plot = NEP_pH_plot, filename = here("Output", "NEP_Plots", "NEP_pH.png"), width = 12, height = 9)

NEP_pH_facetplot <- Data1 %>%
  ggplot(aes(x = NEP, y = pH)) +
  geom_point(aes(color = TREATMENT)) + 
  geom_smooth(method = "lm", formula = y~x) +
  labs(x = "NEP", y = "pH") +
  facet_wrap(~ TREATMENT) +
  theme(axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  stat_regline_equation(label.x = 0, label.y = 8.45, size = 5) + 
  stat_cor(label.x = 0, label.y = 8.4, size = 5) +
  scale_color_manual(labels = c("Control", "Algae-Dominated", "Coral-Dominated", 
                                "Rubble/CCA-Dominated"), values = c("blue", "darkgreen", "coral", "tan"))
NEP_pH_facetplot
#ggsave(plot = NEP_pH_facetplot, filename = here("Output", "NEP_Plots", "NEP_pH_faceted.png"), width = 12, height = 9)

NEP_pH_model <- lmer(pH ~ NEP + (1|TANK_NUM), data = Data1)
check_model(NEP_pH_model)
summary(NEP_pH_model)

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



################################################
### Combine biomass data and carb chem data ### 
# read in biomass (AFDW) data 
afdw <- read_csv(here("Data", "Data_Raw", "Growth", "MO24BEAST_AFDW.csv"))

#filter out 'pre' experiment treatments 
afdw_nopre <- afdw %>%
  filter(!TREATMENT == "Pre")

chem_biomass_data <- afdw_nopre %>% # join afdw data frame with chem data 
  right_join(Data) %>%
  select(CORAL_NUM, GENOTYPE, TREATMENT, TANK_NUM, 
         TA, deltaTA, NEC, mean_AFDW, mean_tissue_biomass, DATE, TIME, pH)

chem_biomass_data_clean <- chem_biomass_data %>%
  filter(TIME %in% c("12:00:00","21:00:00")) %>% 
  group_by(TREATMENT, DATE, TIME, TANK_NUM) 

chem_biomass_plotdata <- chem_biomass_data_clean %>%
  group_by(TANK_NUM, TREATMENT, mean_tissue_biomass, CORAL_NUM, GENOTYPE) %>%
  summarize(pH_mean = mean(pH, na.rm = TRUE),
            pH_se = sd(pH, na.rm = TRUE)/sqrt(n()))

# create model for mean pH influence on mean coral tissue biomass
pH_biomass_model_rando <- lmer(mean_tissue_biomass ~ pH_mean + (1|GENOTYPE), data = chem_biomass_plotdata)
check_model(pH_biomass_model_rando)
summary(pH_biomass_model_rando)

pH_biomass_model <- lm(mean_tissue_biomass ~ pH_mean, data = chem_biomass_plotdata)
check_model(pH_biomass_model)
summary(pH_biomass_model) 

# remove influential plot 
chem_biomass_plotdata2 <- chem_biomass_plotdata %>% 
  filter(mean_tissue_biomass < 0.00075)

# new stats with outlier removed 
pH_biomass_model_rando2 <- lmer(mean_tissue_biomass ~ pH_mean + (1|GENOTYPE), data = chem_biomass_plotdata2)
check_model(pH_biomass_model_rando2)
summary(pH_biomass_model_rando2)

biomass_meanpH_plot <- chem_biomass_plotdata2 %>%
  ggplot(aes(x = pH_mean, y = mean_tissue_biomass)) + 
  geom_point(aes(color = TREATMENT)) +
  labs(y = expression(bold("Mean Tissue Biomass" ~ (g ~ mL^-1 ~ cm^-2))), x= "Mean pH") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_smooth(method = "lm", formula = y~x) + 
  scale_color_manual(values = c("Algae_Dom" = "darkgreen", "Control" = "blue", "Coral_Dom" = "coral",
                                "Rubble_Dom" = "tan"))
biomass_meanpH_plot
ggsave(plot = biomass_meanpH_plot, filename = here("Output", "biomass_meanpH_plot.png"), width = 9, height = 7)


# mean tissue biomass and mean pH with random effect of tank number 
biomass_meanpH_tank <- lmer(mean_tissue_biomass ~ pH_mean + (1|TANK_NUM), data = chem_biomass_plotdata2)
check_model(biomass_meanpH_tank)
summary(biomass_meanpH_tank) # no sig effect of mean pH on mean tissue biomass

## TA and biomass data ## 
chem_biomass_plotdata <- chem_biomass_data_clean %>%
  group_by(TANK_NUM, TREATMENT, mean_tissue_biomass, CORAL_NUM, GENOTYPE) %>%
  summarize(pH_mean = mean(pH, na.rm = TRUE),
            pH_se = sd(pH, na.rm = TRUE)/sqrt(n()), 
            TA_mean = mean(TA, na.rm = TRUE), 
            TA_se = sd(TA, na.rm = TRUE)/sqrt(n()))
chem_biomass_plotdata

chem_biomass_plotdata <- chem_biomass_plotdata %>%
  filter(mean_tissue_biomass < 0.00075)

chem_biomass_plotdata$TREATMENT <- factor(chem_biomass_plotdata$TREATMENT, levels = c("Control", "Algae_Dom", "Coral_Dom", "Rubble_Dom"))

# plot mean tissue biomass as a function of mean TA
biomass_TA_plot <- chem_biomass_plotdata %>%
  ggplot(aes(x = TA_mean, y = mean_tissue_biomass)) + 
  geom_point(aes(color = TREATMENT)) +
  labs(y = expression(bold("Mean Tissue Biomass" ~ (g ~ mL^-1 ~ cm^-2))), x= expression(bold("Mean Total Alkalinity" ~ (µmol ~ kg^-1)))) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray")) +
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values = c("Algae_Dom" = "#E31A1C", "Control" = "green4", "Coral_Dom" = "dodgerblue2",
                                "Rubble_Dom" = "#6A3D9A"))
biomass_TA_plot

#ggsave(plot = biomass_TA_plot, filename = here("Output", "biomass_TA_plot.png"), width = 9, height = 7)

# create model of mean tissue biomass as a function of mean TA 
TA_biomass_model <- lmer(mean_tissue_biomass ~ TA_mean + (1|TANK_NUM), data = chem_biomass_plotdata)
check_model(TA_biomass_model)
summary(TA_biomass_model)
################################################
